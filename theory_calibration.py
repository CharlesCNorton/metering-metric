from __future__ import annotations

import argparse
import itertools
import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np

from metering_theory import bridge_matched_source_amplitude, static_source_amplitude_from_activity
from three_plus_one_lensing import (
    build_archival_residual_map_on_member_grid,
    build_einstein_component_payloads,
    build_multicomponent_lensing_map,
    build_smoothed_background_component_specs,
    compare_theory_and_archival_residuals,
    component_principal_axis_deg,
    component_principal_axis_ratio,
    comparison_objective,
    load_member_geometry_component_specs,
    optimize_rigid_alignment,
    suggested_map_extent,
)


@dataclass(frozen=True)
class ClusterCalibrationContext:
    key: str
    label: str
    member_table: Path
    component_specs: list[dict]
    selection_metadata: dict[str, float]
    map_extent_arcsec: float
    archival_models: list[dict]


def median_or_nan(values: Iterable[float]) -> float:
    finite = [float(value) for value in values if math.isfinite(float(value))]
    if not finite:
        return float("nan")
    return float(np.median(np.asarray(finite, dtype=float)))


def resolve_optional_source_kind(value: object) -> str | None:
    if value is None:
        return None
    rendered = str(value).strip().lower()
    if rendered in {"", "inherit", "none"}:
        return None
    return rendered


def resolve_optional_float(mapping: dict[str, object], key: str) -> float | None:
    value = mapping.get(key, None)
    if value is None:
        return None
    return float(value)


def axis_shift_components(
    shift_parallel_arcsec: float,
    shift_perpendicular_arcsec: float,
    principal_axis_deg: float,
) -> tuple[float, float]:
    angle = math.radians(principal_axis_deg)
    cos_angle = math.cos(angle)
    sin_angle = math.sin(angle)
    shift_x = shift_parallel_arcsec * cos_angle - shift_perpendicular_arcsec * sin_angle
    shift_y = shift_parallel_arcsec * sin_angle + shift_perpendicular_arcsec * cos_angle
    return float(shift_x), float(shift_y)


def axis_coordinates(x_arcsec: float, y_arcsec: float, principal_axis_deg: float) -> tuple[float, float]:
    angle = math.radians(principal_axis_deg)
    cos_angle = math.cos(angle)
    sin_angle = math.sin(angle)
    parallel = x_arcsec * cos_angle + y_arcsec * sin_angle
    perpendicular = -x_arcsec * sin_angle + y_arcsec * cos_angle
    return float(parallel), float(perpendicular)


def from_axis_coordinates(parallel_arcsec: float, perpendicular_arcsec: float, principal_axis_deg: float) -> tuple[float, float]:
    angle = math.radians(principal_axis_deg)
    cos_angle = math.cos(angle)
    sin_angle = math.sin(angle)
    x_arcsec = parallel_arcsec * cos_angle - perpendicular_arcsec * sin_angle
    y_arcsec = parallel_arcsec * sin_angle + perpendicular_arcsec * cos_angle
    return float(x_arcsec), float(y_arcsec)


def anisotropically_shrink_component_specs(
    component_specs: list[dict],
    principal_axis_deg: float,
    shrink_parallel: float,
    shrink_perpendicular: float,
) -> list[dict]:
    clamped_parallel = min(max(shrink_parallel, 0.0), 1.0)
    clamped_perpendicular = min(max(shrink_perpendicular, 0.0), 1.0)
    transformed = []
    for spec in component_specs:
        parallel, perpendicular = axis_coordinates(
            x_arcsec=float(spec["center_x"]),
            y_arcsec=float(spec["center_y"]),
            principal_axis_deg=principal_axis_deg,
        )
        new_x, new_y = from_axis_coordinates(
            parallel_arcsec=(1.0 - clamped_parallel) * parallel,
            perpendicular_arcsec=(1.0 - clamped_perpendicular) * perpendicular,
            principal_axis_deg=principal_axis_deg,
        )
        transformed.append(
            {
                **spec,
                "center_x": new_x,
                "center_y": new_y,
            }
        )
    return transformed


def split_axis_offsets(
    principal_axis_deg: float,
    split_parallel_arcsec: float,
    split_perpendicular_arcsec: float,
) -> list[tuple[float, float]]:
    parallel_half = 0.5 * split_parallel_arcsec
    perpendicular_half = 0.5 * split_perpendicular_arcsec
    if abs(parallel_half) <= 1.0e-12 and abs(perpendicular_half) <= 1.0e-12:
        return [(0.0, 0.0)]
    parallels = [0.0] if abs(parallel_half) <= 1.0e-12 else [-parallel_half, parallel_half]
    perpendiculars = [0.0] if abs(perpendicular_half) <= 1.0e-12 else [-perpendicular_half, perpendicular_half]
    offsets = []
    for parallel in parallels:
        for perpendicular in perpendiculars:
            offsets.append(
                axis_shift_components(
                    shift_parallel_arcsec=parallel,
                    shift_perpendicular_arcsec=perpendicular,
                    principal_axis_deg=principal_axis_deg,
                )
            )
    return offsets


def shift_component_specs(
    component_specs: list[dict],
    shift_x_arcsec: float,
    shift_y_arcsec: float,
    layer_override: str | None = None,
) -> list[dict]:
    shifted = []
    for spec in component_specs:
        shifted.append(
            {
                **spec,
                "center_x": float(spec["center_x"]) + shift_x_arcsec,
                "center_y": float(spec["center_y"]) + shift_y_arcsec,
                "layer": layer_override or spec.get("layer", "member"),
            }
        )
    return shifted


def rotate_component_specs(
    component_specs: list[dict],
    rotation_deg: float,
    layer_override: str | None = None,
) -> list[dict]:
    if abs(rotation_deg) <= 1.0e-12:
        return [
            {
                **spec,
                "layer": layer_override or spec.get("layer", "member"),
            }
            for spec in component_specs
        ]
    angle = math.radians(rotation_deg)
    cos_angle = math.cos(angle)
    sin_angle = math.sin(angle)
    rotated = []
    for spec in component_specs:
        x_arcsec = float(spec["center_x"])
        y_arcsec = float(spec["center_y"])
        row = {
            **spec,
            "center_x": (x_arcsec * cos_angle) - (y_arcsec * sin_angle),
            "center_y": (x_arcsec * sin_angle) + (y_arcsec * cos_angle),
            "layer": layer_override or spec.get("layer", "member"),
        }
        if "orientation_deg" in row:
            row["orientation_deg"] = float(row["orientation_deg"]) + rotation_deg
        rotated.append(row)
    return rotated


def build_shifted_background_component_specs(
    component_specs: list[dict],
    principal_axis_deg: float,
    amplitude_fraction: float,
    size_multiplier: float,
    minimum_size_scale: float,
    position_shrink: float,
    position_shrink_parallel: float | None,
    position_shrink_perpendicular: float | None,
    axis_ratio: float,
    angle_offset_deg: float,
    source_kind: str | None,
    source_beta: float | None,
    source_outer_slope: float | None,
    shift_x_arcsec: float,
    shift_y_arcsec: float,
    layer_name: str,
) -> list[dict]:
    isotropic_shrink = position_shrink
    if position_shrink_parallel is not None or position_shrink_perpendicular is not None:
        isotropic_shrink = 0.0
    specs = build_smoothed_background_component_specs(
        component_specs=component_specs,
        amplitude_fraction=amplitude_fraction,
        size_multiplier=size_multiplier,
        minimum_size_scale=minimum_size_scale,
        position_shrink=isotropic_shrink,
    )
    if position_shrink_parallel is not None or position_shrink_perpendicular is not None:
        specs = anisotropically_shrink_component_specs(
            component_specs=specs,
            principal_axis_deg=principal_axis_deg,
            shrink_parallel=position_shrink if position_shrink_parallel is None else position_shrink_parallel,
            shrink_perpendicular=position_shrink if position_shrink_perpendicular is None else position_shrink_perpendicular,
        )
    shifted = []
    for spec in specs:
        row = {
            **spec,
            "name": f"{layer_name}_{spec['name']}",
            "center_x": float(spec["center_x"]) + shift_x_arcsec,
            "center_y": float(spec["center_y"]) + shift_y_arcsec,
            "axis_ratio": float(axis_ratio),
            "orientation_deg": float(principal_axis_deg + angle_offset_deg),
            "layer": layer_name,
        }
        if source_kind is not None:
            row["source_kind"] = source_kind
        if source_beta is not None:
            row["source_beta"] = float(source_beta)
        if source_outer_slope is not None:
            row["source_outer_slope"] = float(source_outer_slope)
        shifted.append(row)
    return shifted


def discover_model_kappa_paths(input_dir: Path, cluster_prefix: str) -> list[tuple[str, Path]]:
    rows = []
    prefix = f"{cluster_prefix}_"
    suffix = "_kappa.fits"
    for path in sorted(input_dir.glob(f"{cluster_prefix}_*_kappa.fits")):
        name = path.name
        if not (name.startswith(prefix) and name.endswith(suffix)):
            continue
        model_family = name[len(prefix) : -len(suffix)]
        rows.append((model_family, path))
    if not rows:
        raise ValueError(f"No model-family kappa maps found in {input_dir}.")
    return rows


def summarize_model_rows(rows: list[dict]) -> dict:
    return {
        "n_models": int(len(rows)),
        "median_objective": median_or_nan(row["objective"] for row in rows),
        "median_pixel_pearson": median_or_nan(row["pixel_pearson"] for row in rows),
        "median_pixel_spearman": median_or_nan(row["pixel_spearman"] for row in rows),
        "median_resolved_sign_agreement": median_or_nan(row["resolved_sign_agreement"] for row in rows),
        "median_positive_jaccard": median_or_nan(row["positive_jaccard"] for row in rows),
        "median_negative_jaccard": median_or_nan(row["negative_jaccard"] for row in rows),
        "median_harmonic_l2": median_or_nan(row["harmonic_l2"] for row in rows),
        "median_member_score_pearson": median_or_nan(row["member_score_pearson"] for row in rows),
        "median_member_score_spearman": median_or_nan(row["member_score_spearman"] for row in rows),
        "median_member_pass_agreement": median_or_nan(row["member_pass_agreement"] for row in rows),
    }


def load_cluster_context(
    key: str,
    label: str,
    member_table: Path,
    kappa_input_dir: Path,
    top_members: int,
    mass_key: str,
    member_flag_key: str,
    ra_key: str,
    dec_key: str,
    amplitude_scaling_exponent: float,
    size_scaling_exponent: float,
    minimum_size_scale: float,
    auto_extent_padding_arcsec: float,
    minimum_extent_arcsec: float,
    extra_extent_padding_arcsec: float,
    grid_size: int,
    residual_mode: str,
    archive_smooth_sigma_px: float,
    archive_radial_bin_arcsec: float,
) -> ClusterCalibrationContext:
    component_specs, selection_metadata = load_member_geometry_component_specs(
        member_table_path=member_table,
        top_members=top_members,
        mass_key=mass_key,
        member_flag_key=member_flag_key,
        ra_key=ra_key,
        dec_key=dec_key,
        amplitude_scaling_exponent=amplitude_scaling_exponent,
        size_scaling_exponent=size_scaling_exponent,
        minimum_size_scale=minimum_size_scale,
    )
    map_extent = suggested_map_extent(
        component_specs=component_specs,
        padding_arcsec=auto_extent_padding_arcsec + extra_extent_padding_arcsec,
        minimum_extent=minimum_extent_arcsec,
    )
    archival_models = []
    for model_family, kappa_path in discover_model_kappa_paths(kappa_input_dir, key):
        archival_map, archival_metadata, archival_pixel_scale = build_archival_residual_map_on_member_grid(
            kappa_path=kappa_path,
            centroid_ra_deg=float(selection_metadata["centroid_ra_deg"]),
            centroid_dec_deg=float(selection_metadata["centroid_dec_deg"]),
            extent=map_extent,
            grid_size=grid_size,
            residual_mode=residual_mode,
            smooth_sigma_px=archive_smooth_sigma_px,
            radial_bin_arcsec=archive_radial_bin_arcsec,
        )
        archival_models.append(
            {
                "model_family": model_family,
                "kappa_path": str(kappa_path),
                "residual_map": archival_map,
                "metadata": archival_metadata,
                "pixel_scale_arcsec": archival_pixel_scale,
            }
        )
    return ClusterCalibrationContext(
        key=key,
        label=label,
        member_table=member_table,
        component_specs=component_specs,
        selection_metadata=selection_metadata,
        map_extent_arcsec=float(map_extent),
        archival_models=archival_models,
    )


def build_theory_component_specs(
    base_member_specs: list[dict],
    member_shift_parallel_arcsec: float,
    member_shift_perpendicular_arcsec: float,
    member_rotation_deg: float,
    member_source_kind: str | None,
    smooth_background_amplitude_fraction: float,
    smooth_background_size_multiplier: float,
    smooth_background_minimum_size_scale: float,
    smooth_background_position_shrink: float,
    smooth_background_position_shrink_parallel: float | None,
    smooth_background_position_shrink_perpendicular: float | None,
    smooth_background_source_kind: str | None,
    smooth_background_source_beta: float | None,
    smooth_background_source_outer_slope: float | None,
    smooth_background_axis_ratio: float | None,
    smooth_background_angle_offset_deg: float,
    smooth_background_shift_parallel_arcsec: float,
    smooth_background_shift_perpendicular_arcsec: float,
    smooth_background_split_parallel_arcsec: float,
    smooth_background_split_perpendicular_arcsec: float,
    los_amplitude_fraction: float,
    los_size_multiplier: float,
    los_minimum_size_scale: float,
    los_position_shrink: float,
    los_position_shrink_parallel: float | None,
    los_position_shrink_perpendicular: float | None,
    los_source_kind: str | None,
    los_source_beta: float | None,
    los_source_outer_slope: float | None,
    los_axis_ratio: float | None,
    los_angle_offset_deg: float,
    los_shift_parallel_arcsec: float,
    los_shift_perpendicular_arcsec: float,
    los_split_parallel_arcsec: float,
    los_split_perpendicular_arcsec: float,
) -> tuple[list[dict], dict[str, object]]:
    principal_axis_deg = component_principal_axis_deg(base_member_specs)
    member_axis_ratio = component_principal_axis_ratio(base_member_specs)
    resolved_smooth_axis_ratio = member_axis_ratio if smooth_background_axis_ratio is None else smooth_background_axis_ratio
    resolved_los_axis_ratio = member_axis_ratio if los_axis_ratio is None else los_axis_ratio
    member_shift_x, member_shift_y = axis_shift_components(
        shift_parallel_arcsec=member_shift_parallel_arcsec,
        shift_perpendicular_arcsec=member_shift_perpendicular_arcsec,
        principal_axis_deg=principal_axis_deg,
    )
    shifted_member_specs = shift_component_specs(
        component_specs=base_member_specs,
        shift_x_arcsec=member_shift_x,
        shift_y_arcsec=member_shift_y,
    )
    shifted_member_specs = rotate_component_specs(
        component_specs=shifted_member_specs,
        rotation_deg=member_rotation_deg,
    )
    if member_source_kind is not None:
        shifted_member_specs = [
            {
                **spec,
                "source_kind": member_source_kind,
            }
            for spec in shifted_member_specs
        ]
    smooth_background_shift_x, smooth_background_shift_y = axis_shift_components(
        shift_parallel_arcsec=smooth_background_shift_parallel_arcsec,
        shift_perpendicular_arcsec=smooth_background_shift_perpendicular_arcsec,
        principal_axis_deg=principal_axis_deg,
    )
    smooth_background_offsets = split_axis_offsets(
        principal_axis_deg=principal_axis_deg,
        split_parallel_arcsec=smooth_background_split_parallel_arcsec,
        split_perpendicular_arcsec=smooth_background_split_perpendicular_arcsec,
    )
    smooth_weight = 1.0 / max(len(smooth_background_offsets), 1)
    smooth_background_specs = []
    for split_x, split_y in smooth_background_offsets:
        smooth_background_specs.extend(
            build_shifted_background_component_specs(
                component_specs=shifted_member_specs,
                principal_axis_deg=principal_axis_deg,
                amplitude_fraction=smooth_weight * smooth_background_amplitude_fraction,
                size_multiplier=smooth_background_size_multiplier,
                minimum_size_scale=smooth_background_minimum_size_scale,
                position_shrink=smooth_background_position_shrink,
                position_shrink_parallel=smooth_background_position_shrink_parallel,
                position_shrink_perpendicular=smooth_background_position_shrink_perpendicular,
                axis_ratio=resolved_smooth_axis_ratio,
                angle_offset_deg=smooth_background_angle_offset_deg,
                source_kind=smooth_background_source_kind,
                source_beta=smooth_background_source_beta,
                source_outer_slope=smooth_background_source_outer_slope,
                shift_x_arcsec=smooth_background_shift_x + split_x,
                shift_y_arcsec=smooth_background_shift_y + split_y,
                layer_name="smooth_background",
            )
        )
    los_shift_x, los_shift_y = axis_shift_components(
        shift_parallel_arcsec=los_shift_parallel_arcsec,
        shift_perpendicular_arcsec=los_shift_perpendicular_arcsec,
        principal_axis_deg=principal_axis_deg,
    )
    los_offsets = split_axis_offsets(
        principal_axis_deg=principal_axis_deg,
        split_parallel_arcsec=los_split_parallel_arcsec,
        split_perpendicular_arcsec=los_split_perpendicular_arcsec,
    )
    los_weight = 1.0 / max(len(los_offsets), 1)
    los_specs = []
    for split_x, split_y in los_offsets:
        los_specs.extend(
            build_shifted_background_component_specs(
                component_specs=shifted_member_specs,
                principal_axis_deg=principal_axis_deg,
                amplitude_fraction=los_weight * los_amplitude_fraction,
                size_multiplier=los_size_multiplier,
                minimum_size_scale=los_minimum_size_scale,
                position_shrink=los_position_shrink,
                position_shrink_parallel=los_position_shrink_parallel,
                position_shrink_perpendicular=los_position_shrink_perpendicular,
                axis_ratio=resolved_los_axis_ratio,
                angle_offset_deg=los_angle_offset_deg,
                source_kind=los_source_kind,
                source_beta=los_source_beta,
                source_outer_slope=los_source_outer_slope,
                shift_x_arcsec=los_shift_x + split_x,
                shift_y_arcsec=los_shift_y + split_y,
                layer_name="line_of_sight",
            )
        )
    return [*shifted_member_specs, *smooth_background_specs, *los_specs], {
        "principal_axis_deg": float(principal_axis_deg),
        "member_axis_ratio": float(member_axis_ratio),
        "member_shift_x_arcsec": float(member_shift_x),
        "member_shift_y_arcsec": float(member_shift_y),
        "member_rotation_deg": float(member_rotation_deg),
        "member_source_kind": member_source_kind or "inherit",
        "smooth_background_shift_x_arcsec": float(smooth_background_shift_x),
        "smooth_background_shift_y_arcsec": float(smooth_background_shift_y),
        "smooth_background_split_parallel_arcsec": float(smooth_background_split_parallel_arcsec),
        "smooth_background_split_perpendicular_arcsec": float(smooth_background_split_perpendicular_arcsec),
        "smooth_background_position_shrink_parallel": float(
            smooth_background_position_shrink
            if smooth_background_position_shrink_parallel is None
            else smooth_background_position_shrink_parallel
        ),
        "smooth_background_position_shrink_perpendicular": float(
            smooth_background_position_shrink
            if smooth_background_position_shrink_perpendicular is None
            else smooth_background_position_shrink_perpendicular
        ),
        "smooth_background_axis_ratio": float(resolved_smooth_axis_ratio),
        "smooth_background_angle_offset_deg": float(smooth_background_angle_offset_deg),
        "smooth_background_source_kind": smooth_background_source_kind or "inherit",
        "smooth_background_source_beta": None if smooth_background_source_beta is None else float(smooth_background_source_beta),
        "smooth_background_source_outer_slope": None if smooth_background_source_outer_slope is None else float(smooth_background_source_outer_slope),
        "line_of_sight_shift_x_arcsec": float(los_shift_x),
        "line_of_sight_shift_y_arcsec": float(los_shift_y),
        "line_of_sight_split_parallel_arcsec": float(los_split_parallel_arcsec),
        "line_of_sight_split_perpendicular_arcsec": float(los_split_perpendicular_arcsec),
        "line_of_sight_position_shrink_parallel": float(
            los_position_shrink if los_position_shrink_parallel is None else los_position_shrink_parallel
        ),
        "line_of_sight_position_shrink_perpendicular": float(
            los_position_shrink if los_position_shrink_perpendicular is None else los_position_shrink_perpendicular
        ),
        "line_of_sight_axis_ratio": float(resolved_los_axis_ratio),
        "line_of_sight_angle_offset_deg": float(los_angle_offset_deg),
        "line_of_sight_source_kind": los_source_kind or "inherit",
        "line_of_sight_source_beta": None if los_source_beta is None else float(los_source_beta),
        "line_of_sight_source_outer_slope": None if los_source_outer_slope is None else float(los_source_outer_slope),
        "member_component_count": int(len(shifted_member_specs)),
        "smooth_background_component_count": int(len(smooth_background_specs)),
        "line_of_sight_component_count": int(len(los_specs)),
    }


def evaluate_cluster_configuration(
    context: ClusterCalibrationContext,
    source_kind: str,
    source_amplitude: float,
    source_sigma: float,
    source_core_radius: float,
    source_scale_radius: float,
    source_beta: float,
    source_outer_slope: float,
    screening_mass: float,
    density_scale: float,
    gravity_scale: float,
    exterior_match_tolerance: float,
    profile_r_max: float,
    profile_samples: int,
    impact_min: float,
    impact_max: float,
    radial_samples: int,
    r_max: float,
    radial_bin_width: float,
    smooth_sigma_px: float,
    axis_samples: int,
    grid_size: int,
    inner_radius_arcsec: float,
    ring_inner_arcsec: float,
    ring_outer_arcsec: float,
    alignment_max_shift_arcsec: float,
    alignment_shift_steps: int,
    alignment_max_rotation_deg: float,
    alignment_rotation_steps: int,
    component_specs: list[dict],
) -> dict:
    component_payloads = build_einstein_component_payloads(
        component_specs=component_specs,
        source_kind=source_kind,
        source_amplitude=source_amplitude,
        source_sigma=source_sigma,
        source_core_radius=source_core_radius,
        source_scale_radius=source_scale_radius,
        source_beta=source_beta,
        source_outer_slope=source_outer_slope,
        screening_mass=screening_mass,
        density_scale=density_scale,
        gravity_scale=gravity_scale,
        exterior_match_tolerance=exterior_match_tolerance,
        profile_r_max=profile_r_max,
        profile_samples=profile_samples,
        impact_min=impact_min,
        impact_max=impact_max,
        radial_samples=radial_samples,
        r_max=r_max,
    )
    theory_map_payload = build_multicomponent_lensing_map(
        component_payloads=component_payloads,
        extent=context.map_extent_arcsec,
        grid_size=grid_size,
        radial_bin_width=radial_bin_width,
        smooth_sigma_px=smooth_sigma_px,
        axis_samples=axis_samples,
        include_arrays=True,
    )
    theory_map = np.asarray(theory_map_payload["residual_map"], dtype=float)
    member_specs = [spec for spec in component_specs if spec.get("layer", "member") == "member"]
    rows = []
    aligned_rows = []
    for archival_model in context.archival_models:
        comparison = compare_theory_and_archival_residuals(
            theory_map=theory_map,
            archival_map=np.asarray(archival_model["residual_map"], dtype=float),
            extent=context.map_extent_arcsec,
            member_specs=member_specs,
            inner_radius_arcsec=inner_radius_arcsec,
            ring_inner_arcsec=ring_inner_arcsec,
            ring_outer_arcsec=ring_outer_arcsec,
        )
        objective = comparison_objective(comparison)
        rows.append(
            {
                "cluster_key": context.key,
                "model_family": archival_model["model_family"],
                "objective": float(objective),
                "pixel_pearson": float(comparison["pearson_correlation"]),
                "pixel_spearman": float(comparison["spearman_correlation"]),
                "resolved_sign_agreement": float(comparison["resolved_sign_agreement_fraction"]),
                "positive_jaccard": float(comparison["positive_jaccard_resolved"]),
                "negative_jaccard": float(comparison["negative_jaccard_resolved"]),
                "harmonic_l2": float(comparison["harmonic_l2_distance"]),
                "member_score_pearson": float(comparison["member_score_comparison"]["pearson_correlation"]),
                "member_score_spearman": float(comparison["member_score_comparison"]["spearman_correlation"]),
                "member_pass_agreement": float(comparison["member_score_comparison"]["sign_flip_pass_agreement_fraction"]),
            }
        )
        if alignment_max_shift_arcsec > 0.0 or alignment_max_rotation_deg > 0.0:
            alignment = optimize_rigid_alignment(
                theory_map=theory_map,
                archival_map=np.asarray(archival_model["residual_map"], dtype=float),
                extent=context.map_extent_arcsec,
                member_specs=member_specs,
                inner_radius_arcsec=inner_radius_arcsec,
                ring_inner_arcsec=ring_inner_arcsec,
                ring_outer_arcsec=ring_outer_arcsec,
                max_shift_arcsec=alignment_max_shift_arcsec,
                shift_steps=alignment_shift_steps,
                max_rotation_deg=alignment_max_rotation_deg,
                rotation_steps=alignment_rotation_steps,
            )
            aligned = alignment["best_comparison"]
            aligned_rows.append(
                {
                    "cluster_key": context.key,
                    "model_family": archival_model["model_family"],
                    "objective": float(alignment["best_objective"]),
                    "pixel_pearson": float(aligned.get("pearson_correlation", float("nan"))),
                    "pixel_spearman": float(aligned.get("spearman_correlation", float("nan"))),
                    "resolved_sign_agreement": float(aligned.get("resolved_sign_agreement_fraction", float("nan"))),
                    "positive_jaccard": float(aligned.get("positive_jaccard_resolved", float("nan"))),
                    "negative_jaccard": float(aligned.get("negative_jaccard_resolved", float("nan"))),
                    "harmonic_l2": float(aligned.get("harmonic_l2_distance", float("nan"))),
                    "member_score_pearson": float(aligned.get("member_score_comparison", {}).get("pearson_correlation", float("nan"))),
                    "member_score_spearman": float(aligned.get("member_score_comparison", {}).get("spearman_correlation", float("nan"))),
                    "member_pass_agreement": float(
                        aligned.get("member_score_comparison", {}).get("sign_flip_pass_agreement_fraction", float("nan"))
                    ),
                    "shift_x_arcsec": float(alignment["best_shift_x_arcsec"]),
                    "shift_y_arcsec": float(alignment["best_shift_y_arcsec"]),
                    "rotation_deg": float(alignment["best_rotation_deg"]),
                }
            )
    return {
        "cluster_key": context.key,
        "cluster_label": context.label,
        "map_extent_arcsec": float(context.map_extent_arcsec),
        "theory_component_count": int(len(component_specs)),
        "rows": rows,
        "summary": summarize_model_rows(rows),
        "aligned_rows": aligned_rows,
        "aligned_summary": summarize_model_rows(aligned_rows) if aligned_rows else {},
        "theory_map_metrics": {
            "kappa_min": float(theory_map_payload["kappa_min"]),
            "kappa_max": float(theory_map_payload["kappa_max"]),
            "negative_kappa_fraction_resolved": float(theory_map_payload["negative_kappa_fraction_resolved"]),
            "negative_residual_fraction_resolved": float(theory_map_payload["negative_residual_fraction_resolved"]),
            "positive_residual_fraction_resolved": float(theory_map_payload["positive_residual_fraction_resolved"]),
        },
    }


def summarize_global_clusters(cluster_results: list[dict], aligned: bool) -> dict:
    rows = []
    for cluster_result in cluster_results:
        rows.extend(cluster_result["aligned_rows"] if aligned else cluster_result["rows"])
    return summarize_model_rows(rows)


def parameter_product(stage: str, space: dict[str, list[object]]) -> list[dict[str, object]]:
    keys = list(space.keys())
    values = [space[key] for key in keys]
    rows = []
    for combo in itertools.product(*values):
        row = {"stage": stage}
        for key, value in zip(keys, combo):
            row[key] = value
        rows.append(row)
    return rows


def execute_search(
    contexts: list[ClusterCalibrationContext],
    fixed_model_args: dict[str, float | int | str],
    stage_a_space: dict[str, list[object]],
    stage_b_space: dict[str, list[object]],
    top_k_stage_a: int,
) -> dict:
    stage_a_rows = []
    for combo in parameter_product("stage_a", stage_a_space):
        cluster_results = []
        for context in contexts:
            component_specs, structure_metadata = build_theory_component_specs(
                base_member_specs=context.component_specs,
                member_shift_parallel_arcsec=0.0,
                member_shift_perpendicular_arcsec=0.0,
                member_rotation_deg=0.0,
                member_source_kind=resolve_optional_source_kind(combo.get("member_source_kind", None)),
                smooth_background_amplitude_fraction=float(combo["smooth_background_amplitude_fraction"]),
                smooth_background_size_multiplier=float(combo["smooth_background_size_multiplier"]),
                smooth_background_minimum_size_scale=float(combo["smooth_background_minimum_size_scale"]),
                smooth_background_position_shrink=float(combo["smooth_background_position_shrink"]),
                smooth_background_position_shrink_parallel=None,
                smooth_background_position_shrink_perpendicular=None,
                smooth_background_source_kind=resolve_optional_source_kind(combo.get("smooth_background_source_kind", None)),
                smooth_background_source_beta=resolve_optional_float(combo, "smooth_background_source_beta"),
                smooth_background_source_outer_slope=resolve_optional_float(combo, "smooth_background_source_outer_slope"),
                smooth_background_axis_ratio=float(combo.get("smooth_background_axis_ratio", 1.0)),
                smooth_background_angle_offset_deg=float(combo.get("smooth_background_angle_offset_deg", 0.0)),
                smooth_background_shift_parallel_arcsec=0.0,
                smooth_background_shift_perpendicular_arcsec=0.0,
                smooth_background_split_parallel_arcsec=0.0,
                smooth_background_split_perpendicular_arcsec=0.0,
                los_amplitude_fraction=0.0,
                los_size_multiplier=float(combo["los_size_multiplier"]),
                los_minimum_size_scale=float(combo["los_minimum_size_scale"]),
                los_position_shrink=float(combo["los_position_shrink"]),
                los_position_shrink_parallel=None,
                los_position_shrink_perpendicular=None,
                los_source_kind=None,
                los_source_beta=None,
                los_source_outer_slope=None,
                los_axis_ratio=1.0,
                los_angle_offset_deg=0.0,
                los_shift_parallel_arcsec=0.0,
                los_shift_perpendicular_arcsec=0.0,
                los_split_parallel_arcsec=0.0,
                los_split_perpendicular_arcsec=0.0,
            )
            cluster_results.append(
                {
                    **evaluate_cluster_configuration(
                        context=context,
                        source_kind=str(combo["source_kind"]),
                        source_amplitude=float(fixed_model_args["source_amplitude"]),
                        source_sigma=float(fixed_model_args["source_sigma"]),
                        source_core_radius=float(fixed_model_args["source_core_radius"]),
                        source_scale_radius=float(fixed_model_args["source_scale_radius"]),
                        source_beta=float(fixed_model_args["source_beta"]),
                        source_outer_slope=float(fixed_model_args["source_outer_slope"]),
                        screening_mass=float(fixed_model_args["screening_mass"]),
                        density_scale=float(fixed_model_args["density_scale"]),
                        gravity_scale=float(fixed_model_args["gravity_scale"]),
                        exterior_match_tolerance=float(fixed_model_args["exterior_match_tolerance"]),
                        profile_r_max=float(fixed_model_args["profile_r_max"]),
                        profile_samples=int(fixed_model_args["profile_samples"]),
                        impact_min=float(fixed_model_args["impact_min"]),
                        impact_max=float(fixed_model_args["impact_max"]),
                        radial_samples=int(fixed_model_args["radial_samples"]),
                        r_max=float(fixed_model_args["r_max"]),
                        radial_bin_width=float(fixed_model_args["radial_bin_width"]),
                        smooth_sigma_px=float(fixed_model_args["smooth_sigma_px"]),
                        axis_samples=int(fixed_model_args["axis_samples"]),
                        grid_size=int(fixed_model_args["grid_size"]),
                        inner_radius_arcsec=float(fixed_model_args["inner_radius_arcsec"]),
                        ring_inner_arcsec=float(fixed_model_args["ring_inner_arcsec"]),
                        ring_outer_arcsec=float(fixed_model_args["ring_outer_arcsec"]),
                        alignment_max_shift_arcsec=0.0,
                        alignment_shift_steps=1,
                        alignment_max_rotation_deg=0.0,
                        alignment_rotation_steps=1,
                        component_specs=component_specs,
                    ),
                    "structure_metadata": structure_metadata,
                }
            )
        stage_a_rows.append(
            {
                "stage": "stage_a",
                "parameters": combo,
                "cluster_summaries": {
                    cluster_result["cluster_key"]: {
                        **cluster_result["summary"],
                        "structure_metadata": cluster_result["structure_metadata"],
                    }
                    for cluster_result in cluster_results
                },
                "global_summary": summarize_global_clusters(cluster_results, aligned=False),
            }
        )
    stage_a_rows.sort(key=lambda row: row["global_summary"]["median_objective"], reverse=True)

    stage_b_rows = []
    for seed_row in stage_a_rows[:top_k_stage_a]:
        base = seed_row["parameters"]
        for combo in parameter_product("stage_b", stage_b_space):
            parameters = {**base, **combo}
            cluster_results = []
            for context in contexts:
                component_specs, structure_metadata = build_theory_component_specs(
                    base_member_specs=context.component_specs,
                    member_shift_parallel_arcsec=float(parameters["member_shift_parallel_arcsec"]),
                    member_shift_perpendicular_arcsec=float(parameters["member_shift_perpendicular_arcsec"]),
                    member_rotation_deg=float(parameters.get("member_rotation_deg", 0.0)),
                    member_source_kind=resolve_optional_source_kind(parameters.get("member_source_kind", None)),
                    smooth_background_amplitude_fraction=float(parameters["smooth_background_amplitude_fraction"]),
                    smooth_background_size_multiplier=float(parameters["smooth_background_size_multiplier"]),
                    smooth_background_minimum_size_scale=float(parameters["smooth_background_minimum_size_scale"]),
                    smooth_background_position_shrink=float(parameters["smooth_background_position_shrink"]),
                    smooth_background_position_shrink_parallel=resolve_optional_float(
                        parameters,
                        "smooth_background_position_shrink_parallel",
                    ),
                    smooth_background_position_shrink_perpendicular=resolve_optional_float(
                        parameters,
                        "smooth_background_position_shrink_perpendicular",
                    ),
                    smooth_background_source_kind=resolve_optional_source_kind(
                        parameters.get("smooth_background_source_kind", None)
                    ),
                    smooth_background_source_beta=resolve_optional_float(parameters, "smooth_background_source_beta"),
                    smooth_background_source_outer_slope=resolve_optional_float(
                        parameters,
                        "smooth_background_source_outer_slope",
                    ),
                    smooth_background_axis_ratio=float(parameters.get("smooth_background_axis_ratio", 1.0)),
                    smooth_background_angle_offset_deg=float(parameters.get("smooth_background_angle_offset_deg", 0.0)),
                    smooth_background_shift_parallel_arcsec=float(parameters["smooth_background_shift_parallel_arcsec"]),
                    smooth_background_shift_perpendicular_arcsec=float(parameters["smooth_background_shift_perpendicular_arcsec"]),
                    smooth_background_split_parallel_arcsec=float(parameters["smooth_background_split_parallel_arcsec"]),
                    smooth_background_split_perpendicular_arcsec=float(
                        parameters.get("smooth_background_split_perpendicular_arcsec", 0.0)
                    ),
                    los_amplitude_fraction=float(parameters["los_amplitude_fraction"]),
                    los_size_multiplier=float(parameters["los_size_multiplier"]),
                    los_minimum_size_scale=float(parameters["los_minimum_size_scale"]),
                    los_position_shrink=float(parameters["los_position_shrink"]),
                    los_position_shrink_parallel=resolve_optional_float(
                        parameters,
                        "los_position_shrink_parallel",
                    ),
                    los_position_shrink_perpendicular=resolve_optional_float(
                        parameters,
                        "los_position_shrink_perpendicular",
                    ),
                    los_source_kind=resolve_optional_source_kind(parameters.get("los_source_kind", None)),
                    los_source_beta=resolve_optional_float(parameters, "los_source_beta"),
                    los_source_outer_slope=resolve_optional_float(parameters, "los_source_outer_slope"),
                    los_axis_ratio=float(parameters.get("los_axis_ratio", 1.0)),
                    los_angle_offset_deg=float(parameters.get("los_angle_offset_deg", 0.0)),
                    los_shift_parallel_arcsec=float(parameters["los_shift_parallel_arcsec"]),
                    los_shift_perpendicular_arcsec=float(parameters["los_shift_perpendicular_arcsec"]),
                    los_split_parallel_arcsec=float(parameters.get("los_split_parallel_arcsec", 0.0)),
                    los_split_perpendicular_arcsec=float(parameters.get("los_split_perpendicular_arcsec", 0.0)),
                )
                cluster_results.append(
                    {
                        **evaluate_cluster_configuration(
                            context=context,
                            source_kind=str(parameters["source_kind"]),
                            source_amplitude=float(fixed_model_args["source_amplitude"]),
                            source_sigma=float(fixed_model_args["source_sigma"]),
                            source_core_radius=float(fixed_model_args["source_core_radius"]),
                            source_scale_radius=float(fixed_model_args["source_scale_radius"]),
                            source_beta=float(fixed_model_args["source_beta"]),
                            source_outer_slope=float(fixed_model_args["source_outer_slope"]),
                            screening_mass=float(fixed_model_args["screening_mass"]),
                            density_scale=float(fixed_model_args["density_scale"]),
                            gravity_scale=float(fixed_model_args["gravity_scale"]),
                            exterior_match_tolerance=float(fixed_model_args["exterior_match_tolerance"]),
                            profile_r_max=float(fixed_model_args["profile_r_max"]),
                            profile_samples=int(fixed_model_args["profile_samples"]),
                            impact_min=float(fixed_model_args["impact_min"]),
                            impact_max=float(fixed_model_args["impact_max"]),
                            radial_samples=int(fixed_model_args["radial_samples"]),
                            r_max=float(fixed_model_args["r_max"]),
                            radial_bin_width=float(fixed_model_args["radial_bin_width"]),
                            smooth_sigma_px=float(fixed_model_args["smooth_sigma_px"]),
                            axis_samples=int(fixed_model_args["axis_samples"]),
                            grid_size=int(fixed_model_args["grid_size"]),
                            inner_radius_arcsec=float(fixed_model_args["inner_radius_arcsec"]),
                            ring_inner_arcsec=float(fixed_model_args["ring_inner_arcsec"]),
                            ring_outer_arcsec=float(fixed_model_args["ring_outer_arcsec"]),
                            alignment_max_shift_arcsec=0.0,
                            alignment_shift_steps=1,
                            alignment_max_rotation_deg=0.0,
                            alignment_rotation_steps=1,
                            component_specs=component_specs,
                        ),
                        "structure_metadata": structure_metadata,
                    }
                )
            stage_b_rows.append(
                {
                    "stage": "stage_b",
                    "seed_parameters": base,
                    "parameters": parameters,
                    "cluster_summaries": {
                        cluster_result["cluster_key"]: {
                            **cluster_result["summary"],
                            "structure_metadata": cluster_result["structure_metadata"],
                        }
                        for cluster_result in cluster_results
                    },
                    "global_summary": summarize_global_clusters(cluster_results, aligned=False),
                }
            )
    stage_b_rows.sort(key=lambda row: row["global_summary"]["median_objective"], reverse=True)
    best_row = stage_b_rows[0] if stage_b_rows else stage_a_rows[0]
    return {
        "stage_a_rows": stage_a_rows,
        "stage_b_rows": stage_b_rows,
        "best_parameters": best_row["parameters"],
        "best_stage": best_row["stage"],
    }


def final_evaluation(
    contexts: list[ClusterCalibrationContext],
    fixed_model_args: dict[str, float | int | str],
    best_parameters: dict[str, object],
    alignment_max_shift_arcsec: float,
    alignment_shift_steps: int,
    alignment_max_rotation_deg: float,
    alignment_rotation_steps: int,
) -> dict:
    cluster_results = []
    for context in contexts:
        component_specs, structure_metadata = build_theory_component_specs(
            base_member_specs=context.component_specs,
            member_shift_parallel_arcsec=float(best_parameters["member_shift_parallel_arcsec"]),
            member_shift_perpendicular_arcsec=float(best_parameters["member_shift_perpendicular_arcsec"]),
            member_rotation_deg=float(best_parameters.get("member_rotation_deg", 0.0)),
            member_source_kind=resolve_optional_source_kind(best_parameters.get("member_source_kind", None)),
            smooth_background_amplitude_fraction=float(best_parameters["smooth_background_amplitude_fraction"]),
            smooth_background_size_multiplier=float(best_parameters["smooth_background_size_multiplier"]),
            smooth_background_minimum_size_scale=float(best_parameters["smooth_background_minimum_size_scale"]),
            smooth_background_position_shrink=float(best_parameters["smooth_background_position_shrink"]),
            smooth_background_position_shrink_parallel=resolve_optional_float(
                best_parameters,
                "smooth_background_position_shrink_parallel",
            ),
            smooth_background_position_shrink_perpendicular=resolve_optional_float(
                best_parameters,
                "smooth_background_position_shrink_perpendicular",
            ),
            smooth_background_source_kind=resolve_optional_source_kind(best_parameters.get("smooth_background_source_kind", None)),
            smooth_background_source_beta=resolve_optional_float(best_parameters, "smooth_background_source_beta"),
            smooth_background_source_outer_slope=resolve_optional_float(
                best_parameters,
                "smooth_background_source_outer_slope",
            ),
            smooth_background_axis_ratio=float(best_parameters.get("smooth_background_axis_ratio", 1.0)),
            smooth_background_angle_offset_deg=float(best_parameters.get("smooth_background_angle_offset_deg", 0.0)),
            smooth_background_shift_parallel_arcsec=float(best_parameters["smooth_background_shift_parallel_arcsec"]),
            smooth_background_shift_perpendicular_arcsec=float(best_parameters["smooth_background_shift_perpendicular_arcsec"]),
            smooth_background_split_parallel_arcsec=float(best_parameters["smooth_background_split_parallel_arcsec"]),
            smooth_background_split_perpendicular_arcsec=float(
                best_parameters.get("smooth_background_split_perpendicular_arcsec", 0.0)
            ),
            los_amplitude_fraction=float(best_parameters["los_amplitude_fraction"]),
            los_size_multiplier=float(best_parameters["los_size_multiplier"]),
            los_minimum_size_scale=float(best_parameters["los_minimum_size_scale"]),
            los_position_shrink=float(best_parameters["los_position_shrink"]),
            los_position_shrink_parallel=resolve_optional_float(
                best_parameters,
                "los_position_shrink_parallel",
            ),
            los_position_shrink_perpendicular=resolve_optional_float(
                best_parameters,
                "los_position_shrink_perpendicular",
            ),
            los_source_kind=resolve_optional_source_kind(best_parameters.get("los_source_kind", None)),
            los_source_beta=resolve_optional_float(best_parameters, "los_source_beta"),
            los_source_outer_slope=resolve_optional_float(best_parameters, "los_source_outer_slope"),
            los_axis_ratio=float(best_parameters.get("los_axis_ratio", 1.0)),
            los_angle_offset_deg=float(best_parameters.get("los_angle_offset_deg", 0.0)),
            los_shift_parallel_arcsec=float(best_parameters["los_shift_parallel_arcsec"]),
            los_shift_perpendicular_arcsec=float(best_parameters["los_shift_perpendicular_arcsec"]),
            los_split_parallel_arcsec=float(best_parameters.get("los_split_parallel_arcsec", 0.0)),
            los_split_perpendicular_arcsec=float(best_parameters.get("los_split_perpendicular_arcsec", 0.0)),
        )
        cluster_results.append(
            {
                **evaluate_cluster_configuration(
                    context=context,
                    source_kind=str(best_parameters["source_kind"]),
                    source_amplitude=float(fixed_model_args["source_amplitude"]),
                    source_sigma=float(fixed_model_args["source_sigma"]),
                    source_core_radius=float(fixed_model_args["source_core_radius"]),
                    source_scale_radius=float(fixed_model_args["source_scale_radius"]),
                    source_beta=float(fixed_model_args["source_beta"]),
                    source_outer_slope=float(fixed_model_args["source_outer_slope"]),
                    screening_mass=float(fixed_model_args["screening_mass"]),
                    density_scale=float(fixed_model_args["density_scale"]),
                    gravity_scale=float(fixed_model_args["gravity_scale"]),
                    exterior_match_tolerance=float(fixed_model_args["exterior_match_tolerance"]),
                    profile_r_max=float(fixed_model_args["profile_r_max"]),
                    profile_samples=int(fixed_model_args["profile_samples"]),
                    impact_min=float(fixed_model_args["impact_min"]),
                    impact_max=float(fixed_model_args["impact_max"]),
                    radial_samples=int(fixed_model_args["radial_samples"]),
                    r_max=float(fixed_model_args["r_max"]),
                    radial_bin_width=float(fixed_model_args["radial_bin_width"]),
                    smooth_sigma_px=float(fixed_model_args["smooth_sigma_px"]),
                    axis_samples=int(fixed_model_args["axis_samples"]),
                    grid_size=int(fixed_model_args["grid_size"]),
                    inner_radius_arcsec=float(fixed_model_args["inner_radius_arcsec"]),
                    ring_inner_arcsec=float(fixed_model_args["ring_inner_arcsec"]),
                    ring_outer_arcsec=float(fixed_model_args["ring_outer_arcsec"]),
                    alignment_max_shift_arcsec=alignment_max_shift_arcsec,
                    alignment_shift_steps=alignment_shift_steps,
                    alignment_max_rotation_deg=alignment_max_rotation_deg,
                    alignment_rotation_steps=alignment_rotation_steps,
                    component_specs=component_specs,
                ),
                "structure_metadata": structure_metadata,
                "selection_metadata": context.selection_metadata,
            }
        )
    return {
        "best_parameters": best_parameters,
        "cluster_results": cluster_results,
        "global_summary": summarize_global_clusters(cluster_results, aligned=False),
        "aligned_global_summary": summarize_global_clusters(cluster_results, aligned=True),
    }


def write_json(payload: dict, out_path: Path | None) -> None:
    rendered = json.dumps(payload, indent=2)
    if out_path is None:
        print(rendered)
        return
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(rendered + "\n", encoding="utf-8")
    print(rendered)


def command_calibrate_hff_ensemble(args: argparse.Namespace) -> None:
    effective_source_amplitude = float(args.source_amplitude)
    if args.bridge_conversion is not None and args.target_activity_density is not None:
        effective_source_amplitude = static_source_amplitude_from_activity(
            activity_density=float(args.target_activity_density),
            bridge_conversion=float(args.bridge_conversion),
        )
    elif args.reference_activity_density is not None and args.target_activity_density is not None:
        effective_source_amplitude = bridge_matched_source_amplitude(
            reference_source_amplitude=float(args.source_amplitude),
            activity_density=float(args.target_activity_density),
            reference_activity_density=float(args.reference_activity_density),
        )
    smooth_background_source_beta_values = [None] if not args.smooth_background_source_beta_values else args.smooth_background_source_beta_values
    smooth_background_source_outer_slope_values = (
        [None] if not args.smooth_background_source_outer_slope_values else args.smooth_background_source_outer_slope_values
    )
    los_source_beta_values = [None] if not args.los_source_beta_values else args.los_source_beta_values
    los_source_outer_slope_values = [None] if not args.los_source_outer_slope_values else args.los_source_outer_slope_values
    smooth_background_position_shrink_parallel_values = (
        args.smooth_background_position_shrink_parallel_values
        if args.smooth_background_position_shrink_parallel_values
        else args.smooth_background_position_shrink_values
    )
    smooth_background_position_shrink_perpendicular_values = (
        args.smooth_background_position_shrink_perpendicular_values
        if args.smooth_background_position_shrink_perpendicular_values
        else args.smooth_background_position_shrink_values
    )
    los_position_shrink_parallel_values = (
        args.los_position_shrink_parallel_values
        if args.los_position_shrink_parallel_values
        else args.los_position_shrink_values
    )
    los_position_shrink_perpendicular_values = (
        args.los_position_shrink_perpendicular_values
        if args.los_position_shrink_perpendicular_values
        else args.los_position_shrink_values
    )
    extra_extent_padding_arcsec = max(
        abs(value)
        for value in [
            *args.member_shift_parallel_arcsec_values,
            *args.member_shift_perpendicular_arcsec_values,
            *args.smooth_background_shift_parallel_arcsec_values,
            *args.smooth_background_shift_perpendicular_arcsec_values,
            *args.los_shift_parallel_arcsec_values,
            *args.los_shift_perpendicular_arcsec_values,
            *args.smooth_background_split_parallel_arcsec_values,
            *args.smooth_background_split_perpendicular_arcsec_values,
            *args.los_split_parallel_arcsec_values,
            *args.los_split_perpendicular_arcsec_values,
            0.0,
        ]
    ) + 80.0
    coarse_contexts = [
        load_cluster_context(
            key="abell370",
            label="Abell 370",
            member_table=Path(args.abell370_member_table),
            kappa_input_dir=Path(args.abell370_kappa_dir),
            top_members=args.top_members,
            mass_key=args.mass_key,
            member_flag_key=args.member_flag_key,
            ra_key=args.ra_key,
            dec_key=args.dec_key,
            amplitude_scaling_exponent=args.amplitude_scaling_exponent,
            size_scaling_exponent=args.size_scaling_exponent,
            minimum_size_scale=args.minimum_size_scale,
            auto_extent_padding_arcsec=args.auto_extent_padding_arcsec,
            minimum_extent_arcsec=args.minimum_extent_arcsec,
            extra_extent_padding_arcsec=extra_extent_padding_arcsec,
            grid_size=args.search_grid_size,
            residual_mode=args.residual_mode,
            archive_smooth_sigma_px=args.archive_smooth_sigma_px,
            archive_radial_bin_arcsec=args.archive_radial_bin_arcsec,
        ),
        load_cluster_context(
            key="rxcj2248",
            label="RXC J2248 / Abell S1063",
            member_table=Path(args.rxcj2248_member_table),
            kappa_input_dir=Path(args.rxcj2248_kappa_dir),
            top_members=args.top_members,
            mass_key=args.mass_key,
            member_flag_key=args.member_flag_key,
            ra_key=args.ra_key,
            dec_key=args.dec_key,
            amplitude_scaling_exponent=args.amplitude_scaling_exponent,
            size_scaling_exponent=args.size_scaling_exponent,
            minimum_size_scale=args.minimum_size_scale,
            auto_extent_padding_arcsec=args.auto_extent_padding_arcsec,
            minimum_extent_arcsec=args.minimum_extent_arcsec,
            extra_extent_padding_arcsec=extra_extent_padding_arcsec,
            grid_size=args.search_grid_size,
            residual_mode=args.residual_mode,
            archive_smooth_sigma_px=args.archive_smooth_sigma_px,
            archive_radial_bin_arcsec=args.archive_radial_bin_arcsec,
        ),
    ]
    search_result = execute_search(
        contexts=coarse_contexts,
        fixed_model_args={
            "source_amplitude": effective_source_amplitude,
            "source_sigma": args.source_sigma,
            "source_core_radius": args.source_core_radius,
            "source_scale_radius": args.source_scale_radius,
            "source_beta": args.source_beta,
            "source_outer_slope": args.source_outer_slope,
            "screening_mass": args.screening_mass,
            "density_scale": args.density_scale,
            "gravity_scale": args.gravity_scale,
            "exterior_match_tolerance": args.exterior_match_tolerance,
            "profile_r_max": args.search_profile_r_max,
            "profile_samples": args.search_profile_samples,
            "impact_min": args.impact_min,
            "impact_max": args.impact_max,
            "radial_samples": args.search_radial_samples,
            "r_max": args.r_max,
            "radial_bin_width": args.radial_bin_width,
            "smooth_sigma_px": args.smooth_sigma_px,
            "axis_samples": args.axis_samples,
            "grid_size": args.search_grid_size,
            "inner_radius_arcsec": args.inner_radius_arcsec,
            "ring_inner_arcsec": args.ring_inner_arcsec,
            "ring_outer_arcsec": args.ring_outer_arcsec,
        },
        stage_a_space={
            "source_kind": args.source_kind_values,
            "member_source_kind": ["inherit"],
            "smooth_background_amplitude_fraction": args.smooth_background_amplitude_fraction_values,
            "smooth_background_size_multiplier": args.smooth_background_size_multiplier_values,
            "smooth_background_minimum_size_scale": [args.smooth_background_minimum_size_scale],
            "smooth_background_position_shrink": args.smooth_background_position_shrink_values,
            "smooth_background_source_kind": args.smooth_background_source_kind_values,
            "smooth_background_source_beta": smooth_background_source_beta_values,
            "smooth_background_source_outer_slope": smooth_background_source_outer_slope_values,
            "smooth_background_axis_ratio": args.smooth_background_axis_ratio_values,
            "smooth_background_angle_offset_deg": args.smooth_background_angle_offset_deg_values,
            "los_size_multiplier": [args.los_size_multiplier_values[0]],
            "los_minimum_size_scale": [args.los_minimum_size_scale],
            "los_position_shrink": [args.los_position_shrink_values[0]],
        },
        stage_b_space={
            "member_shift_parallel_arcsec": args.member_shift_parallel_arcsec_values,
            "member_shift_perpendicular_arcsec": args.member_shift_perpendicular_arcsec_values,
            "member_rotation_deg": args.member_rotation_deg_values,
            "smooth_background_shift_parallel_arcsec": args.smooth_background_shift_parallel_arcsec_values,
            "smooth_background_shift_perpendicular_arcsec": args.smooth_background_shift_perpendicular_arcsec_values,
            "smooth_background_split_parallel_arcsec": args.smooth_background_split_parallel_arcsec_values,
            "smooth_background_split_perpendicular_arcsec": args.smooth_background_split_perpendicular_arcsec_values,
            "smooth_background_position_shrink_parallel": smooth_background_position_shrink_parallel_values,
            "smooth_background_position_shrink_perpendicular": smooth_background_position_shrink_perpendicular_values,
            "los_amplitude_fraction": args.los_amplitude_fraction_values,
            "los_size_multiplier": args.los_size_multiplier_values,
            "los_minimum_size_scale": [args.los_minimum_size_scale],
            "los_position_shrink": args.los_position_shrink_values,
            "los_position_shrink_parallel": los_position_shrink_parallel_values,
            "los_position_shrink_perpendicular": los_position_shrink_perpendicular_values,
            "los_source_kind": args.los_source_kind_values,
            "los_source_beta": los_source_beta_values,
            "los_source_outer_slope": los_source_outer_slope_values,
            "los_axis_ratio": args.los_axis_ratio_values,
            "los_angle_offset_deg": args.los_angle_offset_deg_values,
            "los_shift_parallel_arcsec": args.los_shift_parallel_arcsec_values,
            "los_shift_perpendicular_arcsec": args.los_shift_perpendicular_arcsec_values,
            "los_split_parallel_arcsec": args.los_split_parallel_arcsec_values,
            "los_split_perpendicular_arcsec": args.los_split_perpendicular_arcsec_values,
        },
        top_k_stage_a=args.top_k_stage_a,
    )
    final_contexts = [
        load_cluster_context(
            key="abell370",
            label="Abell 370",
            member_table=Path(args.abell370_member_table),
            kappa_input_dir=Path(args.abell370_kappa_dir),
            top_members=args.top_members,
            mass_key=args.mass_key,
            member_flag_key=args.member_flag_key,
            ra_key=args.ra_key,
            dec_key=args.dec_key,
            amplitude_scaling_exponent=args.amplitude_scaling_exponent,
            size_scaling_exponent=args.size_scaling_exponent,
            minimum_size_scale=args.minimum_size_scale,
            auto_extent_padding_arcsec=args.auto_extent_padding_arcsec,
            minimum_extent_arcsec=args.minimum_extent_arcsec,
            extra_extent_padding_arcsec=extra_extent_padding_arcsec,
            grid_size=args.final_grid_size,
            residual_mode=args.residual_mode,
            archive_smooth_sigma_px=args.archive_smooth_sigma_px,
            archive_radial_bin_arcsec=args.archive_radial_bin_arcsec,
        ),
        load_cluster_context(
            key="rxcj2248",
            label="RXC J2248 / Abell S1063",
            member_table=Path(args.rxcj2248_member_table),
            kappa_input_dir=Path(args.rxcj2248_kappa_dir),
            top_members=args.top_members,
            mass_key=args.mass_key,
            member_flag_key=args.member_flag_key,
            ra_key=args.ra_key,
            dec_key=args.dec_key,
            amplitude_scaling_exponent=args.amplitude_scaling_exponent,
            size_scaling_exponent=args.size_scaling_exponent,
            minimum_size_scale=args.minimum_size_scale,
            auto_extent_padding_arcsec=args.auto_extent_padding_arcsec,
            minimum_extent_arcsec=args.minimum_extent_arcsec,
            extra_extent_padding_arcsec=extra_extent_padding_arcsec,
            grid_size=args.final_grid_size,
            residual_mode=args.residual_mode,
            archive_smooth_sigma_px=args.archive_smooth_sigma_px,
            archive_radial_bin_arcsec=args.archive_radial_bin_arcsec,
        ),
    ]
    final_result = final_evaluation(
        contexts=final_contexts,
        fixed_model_args={
            "source_amplitude": effective_source_amplitude,
            "source_sigma": args.source_sigma,
            "source_core_radius": args.source_core_radius,
            "source_scale_radius": args.source_scale_radius,
            "source_beta": args.source_beta,
            "source_outer_slope": args.source_outer_slope,
            "screening_mass": args.screening_mass,
            "density_scale": args.density_scale,
            "gravity_scale": args.gravity_scale,
            "exterior_match_tolerance": args.exterior_match_tolerance,
            "profile_r_max": args.final_profile_r_max,
            "profile_samples": args.final_profile_samples,
            "impact_min": args.impact_min,
            "impact_max": args.impact_max,
            "radial_samples": args.final_radial_samples,
            "r_max": args.r_max,
            "radial_bin_width": args.radial_bin_width,
            "smooth_sigma_px": args.smooth_sigma_px,
            "axis_samples": args.axis_samples,
            "grid_size": args.final_grid_size,
            "inner_radius_arcsec": args.inner_radius_arcsec,
            "ring_inner_arcsec": args.ring_inner_arcsec,
            "ring_outer_arcsec": args.ring_outer_arcsec,
        },
        best_parameters=search_result["best_parameters"],
        alignment_max_shift_arcsec=args.alignment_max_shift_arcsec,
        alignment_shift_steps=args.alignment_shift_steps,
        alignment_max_rotation_deg=args.alignment_max_rotation_deg,
        alignment_rotation_steps=args.alignment_rotation_steps,
    )
    payload = {
        "search_inputs": {
            "abell370_member_table": args.abell370_member_table,
            "abell370_kappa_dir": args.abell370_kappa_dir,
            "rxcj2248_member_table": args.rxcj2248_member_table,
            "rxcj2248_kappa_dir": args.rxcj2248_kappa_dir,
            "top_members": args.top_members,
            "amplitude_scaling_exponent": args.amplitude_scaling_exponent,
            "size_scaling_exponent": args.size_scaling_exponent,
            "source_amplitude": args.source_amplitude,
            "effective_source_amplitude": effective_source_amplitude,
            "bridge_conversion": args.bridge_conversion,
            "reference_activity_density": args.reference_activity_density,
            "target_activity_density": args.target_activity_density,
            "density_scale": args.density_scale,
            "gravity_scale": args.gravity_scale,
            "smooth_background_shift_parallel_arcsec_values": args.smooth_background_shift_parallel_arcsec_values,
            "smooth_background_shift_perpendicular_arcsec_values": args.smooth_background_shift_perpendicular_arcsec_values,
            "smooth_background_split_parallel_arcsec_values": args.smooth_background_split_parallel_arcsec_values,
            "smooth_background_split_perpendicular_arcsec_values": args.smooth_background_split_perpendicular_arcsec_values,
            "smooth_background_position_shrink_parallel_values": smooth_background_position_shrink_parallel_values,
            "smooth_background_position_shrink_perpendicular_values": smooth_background_position_shrink_perpendicular_values,
            "smooth_background_source_kind_values": args.smooth_background_source_kind_values,
            "smooth_background_source_beta_values": smooth_background_source_beta_values,
            "smooth_background_source_outer_slope_values": smooth_background_source_outer_slope_values,
            "smooth_background_axis_ratio_values": args.smooth_background_axis_ratio_values,
            "smooth_background_angle_offset_deg_values": args.smooth_background_angle_offset_deg_values,
            "member_rotation_deg_values": args.member_rotation_deg_values,
            "los_split_parallel_arcsec_values": args.los_split_parallel_arcsec_values,
            "los_split_perpendicular_arcsec_values": args.los_split_perpendicular_arcsec_values,
            "los_position_shrink_parallel_values": los_position_shrink_parallel_values,
            "los_position_shrink_perpendicular_values": los_position_shrink_perpendicular_values,
            "los_source_kind_values": args.los_source_kind_values,
            "los_source_beta_values": los_source_beta_values,
            "los_source_outer_slope_values": los_source_outer_slope_values,
            "los_axis_ratio_values": args.los_axis_ratio_values,
            "los_angle_offset_deg_values": args.los_angle_offset_deg_values,
            "residual_mode": args.residual_mode,
            "search_grid_size": args.search_grid_size,
            "final_grid_size": args.final_grid_size,
        },
        "stage_a_best_rows": search_result["stage_a_rows"][: min(args.top_k_stage_a, len(search_result["stage_a_rows"]))],
        "stage_b_best_rows": search_result["stage_b_rows"][: min(args.top_k_stage_b, len(search_result["stage_b_rows"]))],
        "best_stage": search_result["best_stage"],
        "best_parameters": search_result["best_parameters"],
        "final_result": final_result,
    }
    write_json(payload, Path(args.out) if args.out else None)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Calibrate the 3+1D Einstein branch against the HFF model ensemble.")
    subparsers = parser.add_subparsers(dest="command", required=True)
    calibrate = subparsers.add_parser(
        "calibrate-hff-ensemble",
        help="Run a staged search over member-plus-envelope Einstein-branch models against the local HFF ensemble.",
    )
    calibrate.add_argument(
        "--abell370-member-table",
        default="D:\\metering-metric\\combined_firstpass_bandpass\\abell370\\abell370_firstpass_rows.csv",
    )
    calibrate.add_argument(
        "--abell370-kappa-dir",
        default="D:\\metering-metric\\full_experiment_bandpass\\abell370_model_ensemble\\inputs",
    )
    calibrate.add_argument(
        "--rxcj2248-member-table",
        default="D:\\metering-metric\\combined_firstpass_bandpass\\rxcj2248\\rxcj2248_firstpass_rows.csv",
    )
    calibrate.add_argument(
        "--rxcj2248-kappa-dir",
        default="D:\\metering-metric\\full_experiment_bandpass\\rxcj2248_model_ensemble\\inputs",
    )
    calibrate.add_argument("--top-members", type=int, default=40)
    calibrate.add_argument("--mass-key", default="mass_neb")
    calibrate.add_argument("--member-flag-key", default="is_cluster_member")
    calibrate.add_argument("--ra-key", default="ra")
    calibrate.add_argument("--dec-key", default="dec")
    calibrate.add_argument("--amplitude-scaling-exponent", type=float, default=1.0)
    calibrate.add_argument("--size-scaling-exponent", type=float, default=0.5)
    calibrate.add_argument("--minimum-size-scale", type=float, default=0.2)
    calibrate.add_argument("--source-amplitude", type=float, default=1.0e-3)
    calibrate.add_argument("--bridge-conversion", type=float, default=None)
    calibrate.add_argument("--reference-activity-density", type=float, default=None)
    calibrate.add_argument("--target-activity-density", type=float, default=None)
    calibrate.add_argument("--source-sigma", type=float, default=120.0)
    calibrate.add_argument("--source-core-radius", type=float, default=80.0)
    calibrate.add_argument("--source-scale-radius", type=float, default=250.0)
    calibrate.add_argument("--source-beta", type=float, default=0.8)
    calibrate.add_argument("--source-outer-slope", type=float, default=3.0)
    calibrate.add_argument("--screening-mass", type=float, default=0.01)
    calibrate.add_argument("--density-scale", type=float, default=1.0e-6)
    calibrate.add_argument("--gravity-scale", type=float, default=1.0)
    calibrate.add_argument("--exterior-match-tolerance", type=float, default=1.0e-6)
    calibrate.add_argument("--search-profile-r-max", type=float, default=6000.0)
    calibrate.add_argument("--search-profile-samples", type=int, default=1024)
    calibrate.add_argument("--search-radial-samples", type=int, default=40)
    calibrate.add_argument("--final-profile-r-max", type=float, default=8000.0)
    calibrate.add_argument("--final-profile-samples", type=int, default=4096)
    calibrate.add_argument("--final-radial-samples", type=int, default=64)
    calibrate.add_argument("--impact-min", type=float, default=30.0)
    calibrate.add_argument("--impact-max", type=float, default=4000.0)
    calibrate.add_argument("--r-max", type=float, default=80000.0)
    calibrate.add_argument("--radial-bin-width", type=float, default=120.0)
    calibrate.add_argument("--smooth-sigma-px", type=float, default=3.0)
    calibrate.add_argument("--axis-samples", type=int, default=257)
    calibrate.add_argument("--search-grid-size", type=int, default=129)
    calibrate.add_argument("--final-grid-size", type=int, default=513)
    calibrate.add_argument(
        "--source-kind-values",
        nargs="+",
        default=["softened_nfw", "beta_model"],
        choices=["gaussian", "beta_model", "softened_nfw"],
    )
    calibrate.add_argument("--smooth-background-amplitude-fraction-values", nargs="+", type=float, default=[0.25, 0.35, 0.45])
    calibrate.add_argument("--smooth-background-size-multiplier-values", nargs="+", type=float, default=[4.0, 5.0, 6.0])
    calibrate.add_argument("--smooth-background-position-shrink-values", nargs="+", type=float, default=[0.6, 0.75, 0.9])
    calibrate.add_argument("--smooth-background-position-shrink-parallel-values", nargs="*", type=float, default=[])
    calibrate.add_argument("--smooth-background-position-shrink-perpendicular-values", nargs="*", type=float, default=[])
    calibrate.add_argument(
        "--smooth-background-source-kind-values",
        nargs="+",
        default=["inherit", "beta_model", "softened_nfw"],
        choices=["inherit", "gaussian", "beta_model", "softened_nfw"],
    )
    calibrate.add_argument("--smooth-background-source-beta-values", nargs="*", type=float, default=[])
    calibrate.add_argument("--smooth-background-source-outer-slope-values", nargs="*", type=float, default=[])
    calibrate.add_argument("--smooth-background-axis-ratio-values", nargs="+", type=float, default=[1.0])
    calibrate.add_argument("--smooth-background-angle-offset-deg-values", nargs="+", type=float, default=[0.0])
    calibrate.add_argument("--smooth-background-shift-parallel-arcsec-values", nargs="+", type=float, default=[0.0])
    calibrate.add_argument("--smooth-background-shift-perpendicular-arcsec-values", nargs="+", type=float, default=[0.0])
    calibrate.add_argument("--smooth-background-split-parallel-arcsec-values", nargs="+", type=float, default=[0.0])
    calibrate.add_argument("--smooth-background-split-perpendicular-arcsec-values", nargs="+", type=float, default=[0.0])
    calibrate.add_argument("--smooth-background-minimum-size-scale", type=float, default=1.0)
    calibrate.add_argument("--member-shift-parallel-arcsec-values", nargs="+", type=float, default=[0.0, 10.0, -10.0])
    calibrate.add_argument("--member-shift-perpendicular-arcsec-values", nargs="+", type=float, default=[0.0])
    calibrate.add_argument("--member-rotation-deg-values", nargs="+", type=float, default=[0.0])
    calibrate.add_argument("--los-amplitude-fraction-values", nargs="+", type=float, default=[0.0, 0.1, 0.2])
    calibrate.add_argument("--los-size-multiplier-values", nargs="+", type=float, default=[6.0, 8.0])
    calibrate.add_argument("--los-position-shrink-values", nargs="+", type=float, default=[0.4, 0.6])
    calibrate.add_argument("--los-position-shrink-parallel-values", nargs="*", type=float, default=[])
    calibrate.add_argument("--los-position-shrink-perpendicular-values", nargs="*", type=float, default=[])
    calibrate.add_argument(
        "--los-source-kind-values",
        nargs="+",
        default=["inherit", "beta_model", "softened_nfw"],
        choices=["inherit", "gaussian", "beta_model", "softened_nfw"],
    )
    calibrate.add_argument("--los-source-beta-values", nargs="*", type=float, default=[])
    calibrate.add_argument("--los-source-outer-slope-values", nargs="*", type=float, default=[])
    calibrate.add_argument("--los-axis-ratio-values", nargs="+", type=float, default=[1.0])
    calibrate.add_argument("--los-angle-offset-deg-values", nargs="+", type=float, default=[0.0])
    calibrate.add_argument("--los-minimum-size-scale", type=float, default=2.0)
    calibrate.add_argument("--los-shift-parallel-arcsec-values", nargs="+", type=float, default=[0.0, 40.0, -40.0])
    calibrate.add_argument("--los-shift-perpendicular-arcsec-values", nargs="+", type=float, default=[0.0])
    calibrate.add_argument("--los-split-parallel-arcsec-values", nargs="+", type=float, default=[0.0])
    calibrate.add_argument("--los-split-perpendicular-arcsec-values", nargs="+", type=float, default=[0.0])
    calibrate.add_argument("--auto-extent-padding-arcsec", type=float, default=80.0)
    calibrate.add_argument("--minimum-extent-arcsec", type=float, default=120.0)
    calibrate.add_argument(
        "--residual-mode",
        choices=["gaussian_highpass", "radial_median", "radial_median_bandpass"],
        default="radial_median_bandpass",
    )
    calibrate.add_argument("--archive-smooth-sigma-px", type=float, default=18.0)
    calibrate.add_argument("--archive-radial-bin-arcsec", type=float, default=8.0)
    calibrate.add_argument("--inner-radius-arcsec", type=float, default=2.5)
    calibrate.add_argument("--ring-inner-arcsec", type=float, default=5.0)
    calibrate.add_argument("--ring-outer-arcsec", type=float, default=9.0)
    calibrate.add_argument("--top-k-stage-a", type=int, default=3)
    calibrate.add_argument("--top-k-stage-b", type=int, default=5)
    calibrate.add_argument("--alignment-max-shift-arcsec", type=float, default=20.0)
    calibrate.add_argument("--alignment-shift-steps", type=int, default=5)
    calibrate.add_argument("--alignment-max-rotation-deg", type=float, default=20.0)
    calibrate.add_argument("--alignment-rotation-steps", type=int, default=5)
    calibrate.add_argument("--out", type=str, default=None)
    calibrate.set_defaults(func=command_calibrate_hff_ensemble)
    return parser


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
