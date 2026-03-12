"""
Laboratory-side bound workbench for the metering-metric project.

This module turns the repository's Casimir-scale benchmark claim into an
executed, measurement-anchored null bound on the physical coupling alpha.

The benchmark is explicit:

- anomaly anchor: at alpha = 1, the metering-induced pressure correction is
  taken to be 16 times the standard Casimir pressure
- weak-coupling scaling: the anomalous pressure scales as alpha^q with q = 2
  by default
- measurement anchors: use published precision Casimir-pressure uncertainties

This is a benchmark-model laboratory bound, not a final first-principles
derivation of the physical alpha.
"""

from __future__ import annotations

import argparse
import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np
from metering_theory import PLANCK_TIME


HBAR = 1.054_571_817e-34
C = 299_792_458.0
PI = math.pi
NM = 1.0e-9
MICRON = 1.0e-6
MPA = 1.0e-3
UPA = 1.0e-6
NN = 1.0e-9
PN = 1.0e-12


@dataclass(frozen=True)
class AbsolutePressureMeasurement:
    label: str
    absolute_error_pa: float
    separation_min_m: float
    separation_max_m: float
    source: str
    note: str


@dataclass(frozen=True)
class RelativePressureMeasurement:
    label: str
    relative_error: float
    separation_m: float
    source: str
    note: str


DECCA_2003_PRESSURE = AbsolutePressureMeasurement(
    label="Decca et al. dynamic parallel-plate Casimir pressure",
    absolute_error_pa=0.6 * MPA,
    separation_min_m=0.2 * MICRON,
    separation_max_m=1.2 * MICRON,
    source="Phys. Rev. D 68 (2003) 116003 / arXiv:hep-ph/0310157",
    note="Published absolute pressure uncertainty of about 0.6 mPa over 0.2-1.2 um.",
)

DECCA_2007_RELATIVE = RelativePressureMeasurement(
    label="Decca et al. gold-plate Casimir pressure minimum relative error",
    relative_error=0.0019,
    separation_m=160.0 * NM,
    source="Phys. Rev. D 75 (2007) 077101 / arXiv:hep-ph/0703290",
    note="Published minimum relative experimental error of about 0.19% in the 160-750 nm range.",
)

DECCA_2003_STATIC_FORCE_ERROR_PN = 0.3

DEFAULT_ALPHA_TARGETS = (1.0e-1, 1.0e-2, 1.0e-3, 1.0e-4)
DEFAULT_DIFFERENTIAL_ALPHA_TARGETS = (1.0e-2, 1.0e-3, 1.0e-4)
DEFAULT_DIFFERENTIAL_SEPARATIONS_UM = (0.2, 0.3, 0.5)
DEFAULT_PLANAR_SEPARATIONS_UM = (0.2, 0.3, 0.5)
DEFAULT_PLANAR_ALPHA_TARGETS = (1.0e-2, 1.0e-3, 1.0e-4)
DEFAULT_REFERENCE_LAB_SOURCE_DENSITY = 1.0e40
DEFAULT_SOURCE_BRIDGE_LENGTHS_UM = (0.05, 0.1, 0.2, 0.5, 1.0)
DEFAULT_SOURCE_BRIDGE_TIMES_S = (PLANCK_TIME, 1.0e-21, 1.0e-18, 1.0e-15, 1.0e-12)


def perfect_conductor_casimir_pressure(separation_m: float) -> float:
    if separation_m <= 0.0:
        raise ValueError("separation_m must be positive.")
    return (PI * PI * HBAR * C) / (240.0 * separation_m**4)


def benchmark_anomalous_pressure(
    alpha: float,
    separation_m: float,
    anomaly_factor_at_alpha_one: float,
    scaling_power: float,
) -> float:
    if alpha < 0.0:
        raise ValueError("alpha must be nonnegative.")
    if anomaly_factor_at_alpha_one <= 0.0:
        raise ValueError("anomaly_factor_at_alpha_one must be positive.")
    if scaling_power <= 0.0:
        raise ValueError("scaling_power must be positive.")
    casimir_pressure = perfect_conductor_casimir_pressure(separation_m)
    return anomaly_factor_at_alpha_one * casimir_pressure * alpha**scaling_power


def alpha_upper_bound_from_null(
    absolute_error_pa: float,
    separation_m: float,
    anomaly_factor_at_alpha_one: float,
    scaling_power: float,
) -> float:
    if absolute_error_pa <= 0.0:
        raise ValueError("absolute_error_pa must be positive.")
    reference_pressure = anomaly_factor_at_alpha_one * perfect_conductor_casimir_pressure(separation_m)
    return (absolute_error_pa / reference_pressure) ** (1.0 / scaling_power)


def alpha_upper_bound_from_relative_null(
    relative_error: float,
    anomaly_factor_at_alpha_one: float,
    scaling_power: float,
) -> float:
    if relative_error <= 0.0:
        raise ValueError("relative_error must be positive.")
    if anomaly_factor_at_alpha_one <= 0.0:
        raise ValueError("anomaly_factor_at_alpha_one must be positive.")
    if scaling_power <= 0.0:
        raise ValueError("scaling_power must be positive.")
    return (relative_error / anomaly_factor_at_alpha_one) ** (1.0 / scaling_power)


def build_absolute_scan_report(
    measurement: AbsolutePressureMeasurement,
    anomaly_factor_at_alpha_one: float,
    scaling_power: float,
    samples: int,
) -> dict:
    if samples < 2:
        raise ValueError("samples must be at least 2.")
    separations = np.geomspace(measurement.separation_min_m, measurement.separation_max_m, samples)
    rows = []
    for separation_m in separations:
        casimir_pressure = perfect_conductor_casimir_pressure(float(separation_m))
        alpha_upper = alpha_upper_bound_from_null(
            absolute_error_pa=measurement.absolute_error_pa,
            separation_m=float(separation_m),
            anomaly_factor_at_alpha_one=anomaly_factor_at_alpha_one,
            scaling_power=scaling_power,
        )
        rows.append(
            {
                "separation_m": float(separation_m),
                "separation_um": float(separation_m / MICRON),
                "standard_casimir_pressure_pa": float(casimir_pressure),
                "absolute_measurement_error_pa": float(measurement.absolute_error_pa),
                "relative_error_vs_standard_casimir": float(measurement.absolute_error_pa / casimir_pressure),
                "alpha_upper_bound": float(alpha_upper),
            }
        )

    best_row = min(rows, key=lambda row: row["alpha_upper_bound"])
    weakest_row = max(rows, key=lambda row: row["alpha_upper_bound"])
    reference_alpha = 1.0
    reference_rows = [
        {
            **row,
            "benchmark_anomalous_pressure_at_alpha_one_pa": float(
                benchmark_anomalous_pressure(
                    alpha=reference_alpha,
                    separation_m=row["separation_m"],
                    anomaly_factor_at_alpha_one=anomaly_factor_at_alpha_one,
                    scaling_power=scaling_power,
                )
            ),
        }
        for row in rows
    ]
    return {
        "type": "absolute_pressure_scan",
        "measurement": {
            "label": measurement.label,
            "source": measurement.source,
            "note": measurement.note,
            "absolute_error_pa": float(measurement.absolute_error_pa),
            "separation_min_um": float(measurement.separation_min_m / MICRON),
            "separation_max_um": float(measurement.separation_max_m / MICRON),
        },
        "benchmark_model": {
            "anomaly_factor_at_alpha_one": float(anomaly_factor_at_alpha_one),
            "scaling_power": float(scaling_power),
            "benchmark_statement": (
                "At alpha = 1, the anomalous metering-induced pressure is taken to be "
                f"{anomaly_factor_at_alpha_one:.6g} times the standard Casimir pressure."
            ),
        },
        "best_bound": best_row,
        "weakest_bound": weakest_row,
        "rows": reference_rows,
    }


def build_relative_point_report(
    measurement: RelativePressureMeasurement,
    anomaly_factor_at_alpha_one: float,
    scaling_power: float,
) -> dict:
    standard_pressure = perfect_conductor_casimir_pressure(measurement.separation_m)
    alpha_upper = alpha_upper_bound_from_relative_null(
        relative_error=measurement.relative_error,
        anomaly_factor_at_alpha_one=anomaly_factor_at_alpha_one,
        scaling_power=scaling_power,
    )
    return {
        "type": "relative_pressure_point",
        "measurement": {
            "label": measurement.label,
            "source": measurement.source,
            "note": measurement.note,
            "relative_error": float(measurement.relative_error),
            "separation_m": float(measurement.separation_m),
            "separation_um": float(measurement.separation_m / MICRON),
        },
        "derived": {
            "standard_casimir_pressure_pa": float(standard_pressure),
            "absolute_error_pa": float(measurement.relative_error * standard_pressure),
            "alpha_upper_bound": float(alpha_upper),
            "benchmark_anomalous_pressure_at_alpha_one_pa": float(
                benchmark_anomalous_pressure(
                    alpha=1.0,
                    separation_m=measurement.separation_m,
                    anomaly_factor_at_alpha_one=anomaly_factor_at_alpha_one,
                    scaling_power=scaling_power,
                )
            ),
        },
    }


def build_alpha_target_table(
    reference_separation_m: float,
    anomaly_factor_at_alpha_one: float,
    scaling_power: float,
    alpha_targets: Iterable[float],
    comparison_absolute_error_pa: float,
) -> list[dict]:
    standard_pressure = perfect_conductor_casimir_pressure(reference_separation_m)
    rows = []
    for alpha_target in alpha_targets:
        required_relative = anomaly_factor_at_alpha_one * alpha_target**scaling_power
        required_absolute = required_relative * standard_pressure
        rows.append(
            {
                "alpha_target": float(alpha_target),
                "reference_separation_m": float(reference_separation_m),
                "reference_separation_um": float(reference_separation_m / MICRON),
                "standard_casimir_pressure_pa": float(standard_pressure),
                "required_relative_anomaly": float(required_relative),
                "required_absolute_anomaly_pa": float(required_absolute),
                "required_absolute_anomaly_mpa": float(required_absolute / MPA),
                "decca_2003_error_ratio": float(required_absolute / comparison_absolute_error_pa),
            }
        )
    return rows


def build_differential_design_rows(
    separations_m: Iterable[float],
    alpha_targets: Iterable[float],
    anomaly_factor_at_alpha_one: float,
    scaling_power: float,
    plate_area_m2: float,
    active_fill_fraction: float,
    modulation_contrast: float,
    force_floor_n: float,
) -> list[dict]:
    if plate_area_m2 <= 0.0:
        raise ValueError("plate_area_m2 must be positive.")
    if not (0.0 < active_fill_fraction <= 1.0):
        raise ValueError("active_fill_fraction must lie in (0, 1].")
    if not (0.0 < modulation_contrast <= 1.0):
        raise ValueError("modulation_contrast must lie in (0, 1].")
    if force_floor_n <= 0.0:
        raise ValueError("force_floor_n must be positive.")

    rows = []
    effective_fraction = active_fill_fraction * modulation_contrast
    for separation_m in separations_m:
        standard_pressure = perfect_conductor_casimir_pressure(float(separation_m))
        standard_force = standard_pressure * plate_area_m2
        for alpha_target in alpha_targets:
            signal_fraction = anomaly_factor_at_alpha_one * effective_fraction * alpha_target**scaling_power
            differential_pressure = signal_fraction * standard_pressure
            differential_force = differential_pressure * plate_area_m2
            required_area = force_floor_n / differential_pressure if differential_pressure > 0.0 else math.inf
            allowed_gap_drift = signal_fraction * separation_m / 4.0
            rows.append(
                {
                    "separation_m": float(separation_m),
                    "separation_um": float(separation_m / MICRON),
                    "alpha_target": float(alpha_target),
                    "plate_area_m2": float(plate_area_m2),
                    "plate_area_mm2": float(plate_area_m2 * 1.0e6),
                    "standard_casimir_pressure_pa": float(standard_pressure),
                    "standard_casimir_force_n": float(standard_force),
                    "effective_modulated_fraction": float(effective_fraction),
                    "signal_fraction_vs_standard_casimir": float(signal_fraction),
                    "differential_pressure_pa": float(differential_pressure),
                    "differential_pressure_upa": float(differential_pressure / UPA),
                    "differential_force_n": float(differential_force),
                    "differential_force_pn": float(differential_force / PN),
                    "required_common_mode_rejection": float(1.0 / signal_fraction),
                    "uncompensated_gap_drift_limit_m": float(allowed_gap_drift),
                    "uncompensated_gap_drift_pm": float(allowed_gap_drift / 1.0e-12),
                    "force_floor_n": float(force_floor_n),
                    "force_floor_pn": float(force_floor_n / PN),
                    "minimum_plate_area_for_force_floor_m2": float(required_area),
                    "minimum_plate_area_for_force_floor_mm2": float(required_area * 1.0e6),
                }
            )
    return rows


def build_differential_design_report(
    separations_m: Iterable[float],
    alpha_targets: Iterable[float],
    anomaly_factor_at_alpha_one: float,
    scaling_power: float,
    plate_area_m2: float,
    active_fill_fraction: float,
    modulation_contrast: float,
    force_floor_n: float,
) -> dict:
    rows = build_differential_design_rows(
        separations_m=separations_m,
        alpha_targets=alpha_targets,
        anomaly_factor_at_alpha_one=anomaly_factor_at_alpha_one,
        scaling_power=scaling_power,
        plate_area_m2=plate_area_m2,
        active_fill_fraction=active_fill_fraction,
        modulation_contrast=modulation_contrast,
        force_floor_n=force_floor_n,
    )
    keyed_rows = {(row["separation_um"], row["alpha_target"]): row for row in rows}
    reference_key = min(keyed_rows, key=lambda key: (abs(key[0] - 0.2), abs(math.log10(key[1]) + 3.0)))
    reference_row = keyed_rows[reference_key]
    return {
        "benchmark_model": {
            "anomaly_factor_at_alpha_one": float(anomaly_factor_at_alpha_one),
            "scaling_power": float(scaling_power),
            "statement": (
                "Differential design targets inherit the same benchmark anomaly model used in the Casimir null report."
            ),
        },
        "design_inputs": {
            "plate_area_m2": float(plate_area_m2),
            "plate_area_mm2": float(plate_area_m2 * 1.0e6),
            "active_fill_fraction": float(active_fill_fraction),
            "modulation_contrast": float(modulation_contrast),
            "effective_modulated_fraction": float(active_fill_fraction * modulation_contrast),
            "force_floor_n": float(force_floor_n),
            "force_floor_pn": float(force_floor_n / PN),
            "separations_um": [float(separation_m / MICRON) for separation_m in separations_m],
            "alpha_targets": [float(alpha_target) for alpha_target in alpha_targets],
        },
        "reference_design_point": reference_row,
        "rows": rows,
        "executive_summary": {
            "recommended_reading": (
                "For the benchmark family now in the repository, an active-vs-inert Casimir differential becomes an "
                "O(10 uPa) problem at alpha ~ 1e-3 near 0.2 um."
            ),
            "reference_statement": (
                f"At {reference_row['separation_um']:.3g} um and alpha = {reference_row['alpha_target']:.3g}, "
                f"the benchmark differential is {reference_row['differential_pressure_upa']:.6g} uPa "
                f"and {reference_row['differential_force_pn']:.6g} pN for a "
                f"{reference_row['plate_area_mm2']:.6g} mm^2 active area."
            ),
            "force_context_statement": (
                f"That reference force is {reference_row['differential_force_pn'] / DECCA_2003_STATIC_FORCE_ERROR_PN:.6g} "
                f"times the 0.3 pN static-force absolute error quoted in the Decca 2003 Casimir measurement."
            ),
        },
    }


def solve_tridiagonal(
    lower: np.ndarray,
    diagonal: np.ndarray,
    upper: np.ndarray,
    rhs: np.ndarray,
) -> np.ndarray:
    n = diagonal.size
    if n == 0:
        return np.zeros(0, dtype=float)
    if n == 1:
        return np.array([rhs[0] / diagonal[0]], dtype=float)

    c_prime = np.zeros(n - 1, dtype=float)
    d_prime = np.zeros(n, dtype=float)
    c_prime[0] = upper[0] / diagonal[0]
    d_prime[0] = rhs[0] / diagonal[0]

    for index in range(1, n - 1):
        denominator = diagonal[index] - lower[index - 1] * c_prime[index - 1]
        c_prime[index] = upper[index] / denominator
        d_prime[index] = (rhs[index] - lower[index - 1] * d_prime[index - 1]) / denominator

    denominator = diagonal[-1] - lower[-1] * c_prime[-1]
    d_prime[-1] = (rhs[-1] - lower[-1] * d_prime[-2]) / denominator

    solution = np.zeros(n, dtype=float)
    solution[-1] = d_prime[-1]
    for index in range(n - 2, -1, -1):
        solution[index] = d_prime[index] - c_prime[index] * solution[index + 1]
    return solution


def solve_screened_poisson_1d(
    positions: np.ndarray,
    source: np.ndarray,
    screening_length_m: float,
) -> np.ndarray:
    if screening_length_m <= 0.0:
        raise ValueError("screening_length_m must be positive.")
    if positions.ndim != 1 or source.ndim != 1 or positions.size != source.size:
        raise ValueError("positions and source must be 1D arrays of equal length.")
    if positions.size < 3:
        raise ValueError("At least three grid points are required.")

    dz = float(positions[1] - positions[0])
    if not np.allclose(np.diff(positions), dz):
        raise ValueError("positions must be uniformly spaced.")

    screening_mass = 1.0 / screening_length_m
    interior_size = positions.size - 2
    lower = np.full(interior_size - 1, -1.0 / (dz * dz), dtype=float)
    diagonal = np.full(interior_size, 2.0 / (dz * dz) + screening_mass * screening_mass, dtype=float)
    upper = np.full(interior_size - 1, -1.0 / (dz * dz), dtype=float)
    rhs = source[1:-1].astype(float)

    mu = np.zeros_like(source, dtype=float)
    mu[1:-1] = solve_tridiagonal(lower=lower, diagonal=diagonal, upper=upper, rhs=rhs)
    return mu


def build_double_slab_source(
    positions: np.ndarray,
    gap_m: float,
    plate_thickness_m: float,
    source_strength: float,
) -> np.ndarray:
    if gap_m <= 0.0:
        raise ValueError("gap_m must be positive.")
    if plate_thickness_m <= 0.0:
        raise ValueError("plate_thickness_m must be positive.")
    left_inner = -0.5 * gap_m - plate_thickness_m
    left_outer = -0.5 * gap_m
    right_inner = 0.5 * gap_m
    right_outer = 0.5 * gap_m + plate_thickness_m
    source = np.zeros_like(positions, dtype=float)
    left_mask = (positions >= left_inner) & (positions <= left_outer)
    right_mask = (positions >= right_inner) & (positions <= right_outer)
    source[left_mask | right_mask] = source_strength
    return source


def planar_scalar_normal_pressure(
    mu: np.ndarray,
    positions: np.ndarray,
    screening_length_m: float,
) -> np.ndarray:
    screening_mass = 1.0 / screening_length_m
    mu_prime = np.gradient(mu, positions)
    return 0.5 * (mu_prime * mu_prime - (screening_mass * screening_mass) * (mu * mu))


def planar_pressure_difference(
    positions: np.ndarray,
    normal_pressure: np.ndarray,
    gap_m: float,
    plate_thickness_m: float,
) -> dict:
    gap_mask = np.abs(positions) <= 0.5 * gap_m
    outer_extent = float(np.max(np.abs(positions)))
    exterior_threshold = min(
        outer_extent * 0.85,
        0.5 * gap_m + plate_thickness_m + 0.15 * outer_extent,
    )
    exterior_mask = np.abs(positions) >= exterior_threshold
    if np.count_nonzero(exterior_mask) < 8:
        exterior_mask = np.abs(positions) >= 0.8 * outer_extent

    gap_pressure = float(np.mean(normal_pressure[gap_mask]))
    exterior_pressure = float(np.mean(normal_pressure[exterior_mask]))
    signed_difference = exterior_pressure - gap_pressure
    return {
        "gap_pressure": gap_pressure,
        "exterior_pressure": exterior_pressure,
        "signed_pressure_difference": signed_difference,
        "pressure_magnitude": abs(signed_difference),
        "pressure_sign": 0 if signed_difference == 0.0 else int(math.copysign(1.0, signed_difference)),
    }


def planar_double_slab_pressure_base_units(
    gap_m: float,
    plate_thickness_m: float,
    screening_length_m: float,
    source_strength: float,
) -> float:
    if gap_m <= 0.0:
        raise ValueError("gap_m must be positive.")
    if plate_thickness_m <= 0.0:
        raise ValueError("plate_thickness_m must be positive.")
    if screening_length_m <= 0.0:
        raise ValueError("screening_length_m must be positive.")
    screening_mass = 1.0 / screening_length_m
    slab_factor = 1.0 - math.exp(-screening_mass * plate_thickness_m)
    return (
        source_strength
        * source_strength
        * math.exp(-screening_mass * gap_m)
        * slab_factor
        * slab_factor
        / (2.0 * screening_mass * screening_mass)
    )


def build_planar_slab_geometry_report(
    separations_m: Iterable[float],
    alpha_targets: Iterable[float],
    anomaly_factor_at_alpha_one: float,
    scaling_power: float,
    plate_area_m2: float,
    plate_thickness_m: float,
    screening_length_m: float,
    source_strength: float,
    force_floor_n: float,
    samples: int,
    domain_half_width_m: float,
    reference_separation_m: float,
    reference_lab_source_density: float,
) -> dict:
    if samples < 33 or samples % 2 == 0:
        raise ValueError("samples must be an odd integer >= 33.")
    if domain_half_width_m <= 0.0:
        raise ValueError("domain_half_width_m must be positive.")

    separations = sorted({float(reference_separation_m), *[float(value) for value in separations_m]})
    positions = np.linspace(-domain_half_width_m, domain_half_width_m, samples)
    base_rows = []

    for gap_m in separations:
        source = build_double_slab_source(
            positions=positions,
            gap_m=gap_m,
            plate_thickness_m=plate_thickness_m,
            source_strength=source_strength,
        )
        mu = solve_screened_poisson_1d(
            positions=positions,
            source=source,
            screening_length_m=screening_length_m,
        )
        normal_pressure = planar_scalar_normal_pressure(
            mu=mu,
            positions=positions,
            screening_length_m=screening_length_m,
        )
        pressure_summary = planar_pressure_difference(
            positions=positions,
            normal_pressure=normal_pressure,
            gap_m=gap_m,
            plate_thickness_m=plate_thickness_m,
        )
        analytic_base_pressure = planar_double_slab_pressure_base_units(
            gap_m=gap_m,
            plate_thickness_m=plate_thickness_m,
            screening_length_m=screening_length_m,
            source_strength=source_strength,
        )
        base_rows.append(
            {
                "gap_m": float(gap_m),
                "gap_um": float(gap_m / MICRON),
                "source_integral": float(np.trapezoid(source, positions)),
                "max_mu": float(np.max(mu)),
                "midgap_mu": float(np.interp(0.0, positions, mu)),
                "gap_pressure_base_units": float(pressure_summary["gap_pressure"]),
                "exterior_pressure_base_units": float(pressure_summary["exterior_pressure"]),
                "signed_pressure_difference_base_units": float(pressure_summary["signed_pressure_difference"]),
                "pressure_magnitude_base_units": float(pressure_summary["pressure_magnitude"]),
                "analytic_pressure_magnitude_base_units": float(analytic_base_pressure),
                "numeric_to_analytic_ratio": float(
                    pressure_summary["pressure_magnitude"] / analytic_base_pressure
                    if analytic_base_pressure > 0.0
                    else math.nan
                ),
                "pressure_sign": int(pressure_summary["pressure_sign"]),
            }
        )

    reference_row = min(base_rows, key=lambda row: abs(row["gap_m"] - reference_separation_m))
    reference_base_pressure = max(reference_row["pressure_magnitude_base_units"], 1.0e-30)
    reference_target_pressure = anomaly_factor_at_alpha_one * perfect_conductor_casimir_pressure(reference_separation_m)
    calibrated_density_scale = reference_target_pressure / reference_base_pressure
    effective_source_strength = source_strength * math.sqrt(calibrated_density_scale)
    source_conversion_factor = effective_source_strength / reference_lab_source_density
    naive_causal_conversion = 1.0 / (C * screening_length_m)
    effective_persistence_time_s = source_conversion_factor * screening_length_m * screening_length_m

    design_rows = []
    for row in base_rows:
        alpha_one_pressure = calibrated_density_scale * row["pressure_magnitude_base_units"]
        alpha_one_pressure_analytic = calibrated_density_scale * row["analytic_pressure_magnitude_base_units"]
        pure_benchmark_alpha_one = anomaly_factor_at_alpha_one * perfect_conductor_casimir_pressure(row["gap_m"])
        geometry_modifier = alpha_one_pressure / pure_benchmark_alpha_one if pure_benchmark_alpha_one > 0.0 else math.nan
        analytic_geometry_modifier = (
            alpha_one_pressure_analytic / pure_benchmark_alpha_one if pure_benchmark_alpha_one > 0.0 else math.nan
        )
        for alpha_target in alpha_targets:
            predicted_pressure = alpha_one_pressure * alpha_target**scaling_power
            predicted_pressure_analytic = alpha_one_pressure_analytic * alpha_target**scaling_power
            predicted_force = predicted_pressure * plate_area_m2
            predicted_force_analytic = predicted_pressure_analytic * plate_area_m2
            required_area = force_floor_n / predicted_pressure if predicted_pressure > 0.0 else math.inf
            design_rows.append(
                {
                    "gap_m": row["gap_m"],
                    "gap_um": row["gap_um"],
                    "alpha_target": float(alpha_target),
                    "pressure_sign": row["pressure_sign"],
                    "alpha_one_pressure_pa": float(alpha_one_pressure),
                    "alpha_one_pressure_pa_analytic": float(alpha_one_pressure_analytic),
                    "pure_benchmark_alpha_one_pressure_pa": float(pure_benchmark_alpha_one),
                    "geometry_modifier_vs_pure_benchmark": float(geometry_modifier),
                    "geometry_modifier_vs_pure_benchmark_analytic": float(analytic_geometry_modifier),
                    "predicted_differential_pressure_pa": float(predicted_pressure),
                    "predicted_differential_pressure_pa_analytic": float(predicted_pressure_analytic),
                    "predicted_differential_pressure_upa": float(predicted_pressure / UPA),
                    "predicted_differential_pressure_upa_analytic": float(predicted_pressure_analytic / UPA),
                    "predicted_differential_force_n": float(predicted_force),
                    "predicted_differential_force_n_analytic": float(predicted_force_analytic),
                    "predicted_differential_force_pn": float(predicted_force / PN),
                    "predicted_differential_force_pn_analytic": float(predicted_force_analytic / PN),
                    "plate_area_m2": float(plate_area_m2),
                    "plate_area_mm2": float(plate_area_m2 * 1.0e6),
                    "force_floor_n": float(force_floor_n),
                    "force_floor_pn": float(force_floor_n / PN),
                    "minimum_plate_area_for_force_floor_m2": float(required_area),
                    "minimum_plate_area_for_force_floor_mm2": float(required_area * 1.0e6),
                }
            )

    reference_design = min(
        design_rows,
        key=lambda row: (abs(row["gap_um"] - reference_separation_m / MICRON), abs(math.log10(row["alpha_target"]) + 3.0)),
    )
    return {
        "field_model": {
            "equation": "-mu'' + m_mu^2 mu = J(z)",
            "plate_thickness_m": float(plate_thickness_m),
            "plate_thickness_um": float(plate_thickness_m / MICRON),
            "screening_length_m": float(screening_length_m),
            "screening_length_um": float(screening_length_m / MICRON),
            "source_strength": float(source_strength),
            "domain_half_width_m": float(domain_half_width_m),
            "domain_half_width_um": float(domain_half_width_m / MICRON),
            "samples": int(samples),
            "statement": (
                "This report derives a geometry-dependent slab pressure scale from the screened metering field equation "
                "and the canonical scalar normal-stress profile, then calibrates that scale to the repository's "
                "reference anomaly amplitude at the chosen reference gap."
            ),
        },
        "benchmark_calibration": {
            "reference_separation_m": float(reference_separation_m),
            "reference_separation_um": float(reference_separation_m / MICRON),
            "anomaly_factor_at_alpha_one": float(anomaly_factor_at_alpha_one),
            "scaling_power": float(scaling_power),
            "reference_target_pressure_pa": float(reference_target_pressure),
            "reference_base_pressure": float(reference_base_pressure),
            "calibrated_density_scale": float(calibrated_density_scale),
        },
        "source_normalization_gap": {
            "effective_source_strength_required": float(effective_source_strength),
            "reference_lab_source_density": float(reference_lab_source_density),
            "effective_conversion_factor": float(source_conversion_factor),
            "naive_causal_conversion_factor": float(naive_causal_conversion),
            "conversion_vs_naive_causal": float(source_conversion_factor / naive_causal_conversion),
            "effective_persistence_time_s": float(effective_persistence_time_s),
            "effective_persistence_time_in_planck_units": float(effective_persistence_time_s / PLANCK_TIME),
            "statement": (
                "This is the remaining laboratory normalization gap in one number: the effective planar source strength "
                "needed to reproduce the reference pressure, and the implied conversion from the repository's laboratory "
                "decoherence-rate-density scale to the static screened source entering the slab equation."
            ),
        },
        "closed_form_model": {
            "pressure_formula": "P(g) = J0^2 * exp(-m_mu g) * (1 - exp(-m_mu t))^2 / (2 m_mu^2)",
            "statement": (
                "This is the exact planar double-slab interaction pressure for the linear screened field equation, "
                "before benchmark amplitude calibration."
            ),
        },
        "base_gap_rows": base_rows,
        "design_rows": design_rows,
        "reference_design_point": reference_design,
        "executive_summary": {
            "recommended_reading": (
                "This is the first geometry-dependent laboratory signal model in the repository with an exact "
                "closed-form slab interaction law. It keeps the benchmark amplitude fixed at the reference gap "
                "but replaces the pure 1/a^4 extrapolation with a screened-field interaction shape."
            ),
            "reference_statement": (
                f"At {reference_design['gap_um']:.3g} um and alpha = {reference_design['alpha_target']:.3g}, "
                f"the planar slab model predicts {reference_design['predicted_differential_pressure_upa_analytic']:.6g} uPa "
                f"and {reference_design['predicted_differential_force_pn_analytic']:.6g} pN for a "
                f"{reference_design['plate_area_mm2']:.6g} mm^2 active area."
            ),
            "normalization_statement": (
                f"Matching the reference gap requires an effective slab source strength of {effective_source_strength:.6g}, "
                f"which corresponds to a conversion factor of {source_conversion_factor:.6g} relative to the repository's "
                f"laboratory decoherence-rate-density scale {reference_lab_source_density:.6g}. "
                f"That is {source_conversion_factor / naive_causal_conversion:.6g} times the naive causal screening scale "
                f"1/(c l_s). If read as a simple persistence-time bridge kappa = tau/l_s^2, the implied timescale "
                f"is {effective_persistence_time_s:.6g} s = {effective_persistence_time_s / PLANCK_TIME:.6g} t_P."
            ),
            "force_context_statement": (
                f"The reference differential force {reference_design['predicted_differential_force_pn_analytic']:.6g} pN is "
                f"{reference_design['predicted_differential_force_pn_analytic'] / DECCA_2003_STATIC_FORCE_ERROR_PN:.6g} "
                f"times the 0.3 pN static-force absolute error quoted in the Decca 2003 Casimir measurement."
            ),
        },
    }


def build_source_bridge_report(
    anomaly_factor_at_alpha_one: float,
    scaling_power: float,
    plate_area_m2: float,
    plate_thickness_m: float,
    screening_length_m: float,
    source_strength: float,
    force_floor_n: float,
    samples: int,
    domain_half_width_m: float,
    reference_separation_m: float,
    reference_lab_source_density: float,
    alpha_targets: Iterable[float],
    separations_m: Iterable[float],
    candidate_lengths_m: Iterable[float],
    candidate_times_s: Iterable[float],
) -> dict:
    planar_report = build_planar_slab_geometry_report(
        separations_m=separations_m,
        alpha_targets=alpha_targets,
        anomaly_factor_at_alpha_one=anomaly_factor_at_alpha_one,
        scaling_power=scaling_power,
        plate_area_m2=plate_area_m2,
        plate_thickness_m=plate_thickness_m,
        screening_length_m=screening_length_m,
        source_strength=source_strength,
        force_floor_n=force_floor_n,
        samples=samples,
        domain_half_width_m=domain_half_width_m,
        reference_separation_m=reference_separation_m,
        reference_lab_source_density=reference_lab_source_density,
    )
    normalization_gap = planar_report["source_normalization_gap"]
    required_conversion = float(normalization_gap["effective_conversion_factor"])
    candidate_length_rows = []
    for length_m in candidate_lengths_m:
        required_tau = required_conversion * length_m * length_m
        candidate_length_rows.append(
            {
                "transverse_scale_m": float(length_m),
                "transverse_scale_um": float(length_m / MICRON),
                "required_persistence_time_s_for_unit_efficiency": float(required_tau),
                "required_persistence_time_in_planck_units_for_unit_efficiency": float(required_tau / PLANCK_TIME),
            }
        )

    candidate_time_rows = []
    for time_s in candidate_times_s:
        required_efficiency = required_conversion * screening_length_m * screening_length_m / time_s
        candidate_time_rows.append(
            {
                "persistence_time_s": float(time_s),
                "persistence_time_in_planck_units": float(time_s / PLANCK_TIME),
                "required_efficiency_at_screening_length": float(required_efficiency),
            }
        )

    planck_tick_efficiency = required_conversion * screening_length_m * screening_length_m / PLANCK_TIME
    planck_tick_length_for_unit_efficiency = math.sqrt(PLANCK_TIME / required_conversion)
    screening_length_row = min(candidate_length_rows, key=lambda row: abs(row["transverse_scale_m"] - screening_length_m))
    nearest_unit_efficiency_length_row = min(
        candidate_length_rows,
        key=lambda row: abs(math.log10(max(row["required_persistence_time_in_planck_units_for_unit_efficiency"], 1.0e-300))),
    )
    return {
        "bridge_model": {
            "dynamic_activity_formula": "R(x) = sum_i n_i(x) * gamma_D,i(x)",
            "static_source_formula": "J(x) = kappa_J * R(x)",
            "bridge_formula": "kappa_J = eta_J * tau_p / L_perp^2",
            "dynamic_activity_units": "m^-3 s^-1",
            "static_source_units": "m^-5",
            "bridge_units": "s m^-2",
            "statement": (
                "The laboratory source bridge is formulated here as a reduction from the dynamic metering-activity "
                "density R to the static screened source J. The factor kappa_J is interpreted as a persistence-time "
                "factor divided by a transverse coarse-graining area."
            ),
        },
        "planar_reference": {
            "reference_separation_m": float(reference_separation_m),
            "reference_separation_um": float(reference_separation_m / MICRON),
            "screening_length_m": float(screening_length_m),
            "screening_length_um": float(screening_length_m / MICRON),
            "plate_thickness_m": float(plate_thickness_m),
            "plate_thickness_um": float(plate_thickness_m / MICRON),
            "required_bridge_conversion": float(required_conversion),
            "reference_activity_density": float(reference_lab_source_density),
            "reference_effective_source_strength": float(normalization_gap["effective_source_strength_required"]),
        },
        "candidate_length_rows": candidate_length_rows,
        "candidate_time_rows": candidate_time_rows,
        "planck_tick_screening_candidate": {
            "persistence_time_s": float(PLANCK_TIME),
            "transverse_scale_m": float(screening_length_m),
            "transverse_scale_um": float(screening_length_m / MICRON),
            "required_efficiency_eta_J": float(planck_tick_efficiency),
        },
        "unit_efficiency_planck_tick_candidate": {
            "required_transverse_scale_m": float(planck_tick_length_for_unit_efficiency),
            "required_transverse_scale_um": float(planck_tick_length_for_unit_efficiency / MICRON),
            "persistence_time_s": float(PLANCK_TIME),
            "persistence_time_in_planck_units": 1.0,
        },
        "recommended_bridge_family": {
            "formula": "J(x) = eta_J * (t_P / L_perp(x)^2) * R(x)",
            "screening_length_specialization": "J(x) = eta_J * (t_P / l_s(x)^2) * R(x)",
            "required_eta_J_if_L_perp_equals_l_s": float(planck_tick_efficiency),
            "statement": (
                "At the current slab benchmark point, a Planck-tick occupancy bridge is the first non-naive bridge "
                "family that lands in a dimensionless efficiency range below unity rather than at absurd suppression. "
                "With L_perp = l_s = 0.2 um, the required eta_J is about 1.41e-2."
            ),
        },
        "executive_summary": {
            "recommended_reading": (
                "The exact slab law no longer leaves the laboratory bridge unconstrained. It selects a required "
                "conversion kappa_J and therefore a narrow family of viable bridge laws."
            ),
            "bridge_statement": (
                f"Matching the current slab benchmark requires kappa_J = {required_conversion:.6g} s m^-2. "
                f"If J = eta_J * (t_P / l_s^2) * R with l_s = {screening_length_m / MICRON:.6g} um, "
                f"the required eta_J is {planck_tick_efficiency:.6g}."
            ),
            "scale_statement": (
                f"With eta_J = 1 and tau_p = t_P, the corresponding transverse coarse-graining scale is "
                f"{planck_tick_length_for_unit_efficiency / MICRON:.6g} um."
            ),
            "screening_length_statement": (
                f"For unit efficiency at L_perp = {screening_length_row['transverse_scale_um']:.6g} um, the required "
                f"persistence time is {screening_length_row['required_persistence_time_s_for_unit_efficiency']:.6g} s "
                f"= {screening_length_row['required_persistence_time_in_planck_units_for_unit_efficiency']:.6g} t_P."
            ),
            "nearest_length_statement": (
                f"Among the current candidate transverse scales, the one nearest unit Planck-scale efficiency is "
                f"{nearest_unit_efficiency_length_row['transverse_scale_um']:.6g} um."
            ),
        },
        "embedded_planar_report": planar_report,
    }


def build_casimir_benchmark_report(
    anomaly_factor_at_alpha_one: float,
    scaling_power: float,
    samples: int,
    reference_separation_m: float,
    alpha_targets: Iterable[float],
) -> dict:
    absolute_scan = build_absolute_scan_report(
        measurement=DECCA_2003_PRESSURE,
        anomaly_factor_at_alpha_one=anomaly_factor_at_alpha_one,
        scaling_power=scaling_power,
        samples=samples,
    )
    relative_point = build_relative_point_report(
        measurement=DECCA_2007_RELATIVE,
        anomaly_factor_at_alpha_one=anomaly_factor_at_alpha_one,
        scaling_power=scaling_power,
    )
    strongest_anchor = min(
        [
            {
                "label": absolute_scan["measurement"]["label"],
                "source": absolute_scan["measurement"]["source"],
                "alpha_upper_bound": absolute_scan["best_bound"]["alpha_upper_bound"],
                "anchor_type": absolute_scan["type"],
                "reference_separation_um": absolute_scan["best_bound"]["separation_um"],
            },
            {
                "label": relative_point["measurement"]["label"],
                "source": relative_point["measurement"]["source"],
                "alpha_upper_bound": relative_point["derived"]["alpha_upper_bound"],
                "anchor_type": relative_point["type"],
                "reference_separation_um": relative_point["measurement"]["separation_um"],
            },
        ],
        key=lambda row: row["alpha_upper_bound"],
    )
    alpha_target_table = build_alpha_target_table(
        reference_separation_m=reference_separation_m,
        anomaly_factor_at_alpha_one=anomaly_factor_at_alpha_one,
        scaling_power=scaling_power,
        alpha_targets=alpha_targets,
        comparison_absolute_error_pa=DECCA_2003_PRESSURE.absolute_error_pa,
    )
    return {
        "benchmark_model": {
            "anomaly_factor_at_alpha_one": float(anomaly_factor_at_alpha_one),
            "scaling_power": float(scaling_power),
            "statement": (
                "At alpha = 1, the anomalous metering-induced pressure is benchmarked as "
                f"{anomaly_factor_at_alpha_one:.6g} times the standard Casimir pressure, "
                f"with weak-coupling scaling proportional to alpha^{scaling_power:.6g}."
            ),
        },
        "measurement_anchors": {
            "decca_2003_absolute_scan": absolute_scan,
            "decca_2007_relative_point": relative_point,
        },
        "strongest_bound": strongest_anchor,
        "reference_target_scan": {
            "reference_separation_m": float(reference_separation_m),
            "reference_separation_um": float(reference_separation_m / MICRON),
            "alpha_target_table": alpha_target_table,
        },
        "executive_summary": {
            "benchmark_lab_bound_alpha_upper": float(strongest_anchor["alpha_upper_bound"]),
            "recommended_reading": (
                "This is a benchmark-model null bound. It is stronger than the older internal "
                "WEC consistency cutoff, but it still rests on the repository's present "
                "Casimir-scale anomaly model rather than a completed first-principles lab derivation."
            ),
        },
    }


def write_json_report(report: dict, out_path: str | None) -> None:
    rendered = json.dumps(report, indent=2)
    if out_path is None:
        print(rendered)
        return
    destination = Path(out_path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(rendered + "\n", encoding="utf-8")
    print(rendered)


def command_casimir_benchmark_report(args: argparse.Namespace) -> None:
    report = build_casimir_benchmark_report(
        anomaly_factor_at_alpha_one=args.anomaly_factor_at_alpha_one,
        scaling_power=args.scaling_power,
        samples=args.samples,
        reference_separation_m=args.reference_separation_um * MICRON,
        alpha_targets=args.alpha_targets,
    )
    write_json_report(report, args.out)


def command_casimir_differential_design(args: argparse.Namespace) -> None:
    report = build_differential_design_report(
        separations_m=[value * MICRON for value in args.separations_um],
        alpha_targets=args.alpha_targets,
        anomaly_factor_at_alpha_one=args.anomaly_factor_at_alpha_one,
        scaling_power=args.scaling_power,
        plate_area_m2=args.plate_area_mm2 * 1.0e-6,
        active_fill_fraction=args.active_fill_fraction,
        modulation_contrast=args.modulation_contrast,
        force_floor_n=args.force_floor_pn * PN,
    )
    write_json_report(report, args.out)


def command_casimir_planar_slab_report(args: argparse.Namespace) -> None:
    report = build_planar_slab_geometry_report(
        separations_m=[value * MICRON for value in args.separations_um],
        alpha_targets=args.alpha_targets,
        anomaly_factor_at_alpha_one=args.anomaly_factor_at_alpha_one,
        scaling_power=args.scaling_power,
        plate_area_m2=args.plate_area_mm2 * 1.0e-6,
        plate_thickness_m=args.plate_thickness_um * MICRON,
        screening_length_m=args.screening_length_um * MICRON,
        source_strength=args.source_strength,
        force_floor_n=args.force_floor_pn * PN,
        samples=args.samples,
        domain_half_width_m=args.domain_half_width_um * MICRON,
        reference_separation_m=args.reference_separation_um * MICRON,
        reference_lab_source_density=args.reference_lab_source_density,
    )
    write_json_report(report, args.out)


def command_casimir_source_bridge_report(args: argparse.Namespace) -> None:
    report = build_source_bridge_report(
        anomaly_factor_at_alpha_one=args.anomaly_factor_at_alpha_one,
        scaling_power=args.scaling_power,
        plate_area_m2=args.plate_area_mm2 * 1.0e-6,
        plate_thickness_m=args.plate_thickness_um * MICRON,
        screening_length_m=args.screening_length_um * MICRON,
        source_strength=args.source_strength,
        force_floor_n=args.force_floor_pn * PN,
        samples=args.samples,
        domain_half_width_m=args.domain_half_width_um * MICRON,
        reference_separation_m=args.reference_separation_um * MICRON,
        reference_lab_source_density=args.reference_lab_source_density,
        alpha_targets=args.alpha_targets,
        separations_m=[value * MICRON for value in args.separations_um],
        candidate_lengths_m=[value * MICRON for value in args.bridge_lengths_um],
        candidate_times_s=args.bridge_times_s,
    )
    write_json_report(report, args.out)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Laboratory-side benchmark bounds for the metering-metric coupling.",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    decca = subparsers.add_parser(
        "casimir-benchmark-report",
        help="Build a measurement-anchored Casimir benchmark report for alpha.",
    )
    decca.add_argument("--anomaly-factor-at-alpha-one", type=float, default=16.0)
    decca.add_argument("--scaling-power", type=float, default=2.0)
    decca.add_argument("--samples", type=int, default=64)
    decca.add_argument("--reference-separation-um", type=float, default=0.2)
    decca.add_argument("--alpha-targets", type=float, nargs="+", default=list(DEFAULT_ALPHA_TARGETS))
    decca.add_argument("--out", type=str, default=None)
    decca.set_defaults(func=command_casimir_benchmark_report)

    differential = subparsers.add_parser(
        "casimir-differential-design",
        help="Build an active-vs-inert Casimir differential design target report.",
    )
    differential.add_argument("--anomaly-factor-at-alpha-one", type=float, default=16.0)
    differential.add_argument("--scaling-power", type=float, default=2.0)
    differential.add_argument(
        "--separations-um",
        type=float,
        nargs="+",
        default=list(DEFAULT_DIFFERENTIAL_SEPARATIONS_UM),
    )
    differential.add_argument(
        "--alpha-targets",
        type=float,
        nargs="+",
        default=list(DEFAULT_DIFFERENTIAL_ALPHA_TARGETS),
    )
    differential.add_argument("--plate-area-mm2", type=float, default=1.0)
    differential.add_argument("--active-fill-fraction", type=float, default=1.0)
    differential.add_argument("--modulation-contrast", type=float, default=1.0)
    differential.add_argument("--force-floor-pn", type=float, default=10.0)
    differential.add_argument("--out", type=str, default=None)
    differential.set_defaults(func=command_casimir_differential_design)

    planar = subparsers.add_parser(
        "casimir-planar-slab-report",
        help="Build a geometry-dependent slab-field Casimir signal report.",
    )
    planar.add_argument("--anomaly-factor-at-alpha-one", type=float, default=16.0)
    planar.add_argument("--scaling-power", type=float, default=2.0)
    planar.add_argument("--reference-separation-um", type=float, default=0.2)
    planar.add_argument(
        "--separations-um",
        type=float,
        nargs="+",
        default=list(DEFAULT_PLANAR_SEPARATIONS_UM),
    )
    planar.add_argument(
        "--alpha-targets",
        type=float,
        nargs="+",
        default=list(DEFAULT_PLANAR_ALPHA_TARGETS),
    )
    planar.add_argument("--plate-area-mm2", type=float, default=1.0)
    planar.add_argument("--plate-thickness-um", type=float, default=0.05)
    planar.add_argument("--screening-length-um", type=float, default=0.2)
    planar.add_argument("--source-strength", type=float, default=1.0)
    planar.add_argument("--force-floor-pn", type=float, default=10.0)
    planar.add_argument("--samples", type=int, default=4001)
    planar.add_argument("--domain-half-width-um", type=float, default=2.0)
    planar.add_argument("--reference-lab-source-density", type=float, default=DEFAULT_REFERENCE_LAB_SOURCE_DENSITY)
    planar.add_argument("--out", type=str, default=None)
    planar.set_defaults(func=command_casimir_planar_slab_report)

    bridge = subparsers.add_parser(
        "casimir-source-bridge-report",
        help="Build a source-normalization bridge report on top of the exact slab law.",
    )
    bridge.add_argument("--anomaly-factor-at-alpha-one", type=float, default=16.0)
    bridge.add_argument("--scaling-power", type=float, default=2.0)
    bridge.add_argument("--reference-separation-um", type=float, default=0.2)
    bridge.add_argument(
        "--separations-um",
        type=float,
        nargs="+",
        default=list(DEFAULT_PLANAR_SEPARATIONS_UM),
    )
    bridge.add_argument(
        "--alpha-targets",
        type=float,
        nargs="+",
        default=list(DEFAULT_PLANAR_ALPHA_TARGETS),
    )
    bridge.add_argument("--plate-area-mm2", type=float, default=1.0)
    bridge.add_argument("--plate-thickness-um", type=float, default=0.05)
    bridge.add_argument("--screening-length-um", type=float, default=0.2)
    bridge.add_argument("--source-strength", type=float, default=1.0)
    bridge.add_argument("--force-floor-pn", type=float, default=10.0)
    bridge.add_argument("--samples", type=int, default=4001)
    bridge.add_argument("--domain-half-width-um", type=float, default=2.0)
    bridge.add_argument("--reference-lab-source-density", type=float, default=DEFAULT_REFERENCE_LAB_SOURCE_DENSITY)
    bridge.add_argument(
        "--bridge-lengths-um",
        type=float,
        nargs="+",
        default=list(DEFAULT_SOURCE_BRIDGE_LENGTHS_UM),
    )
    bridge.add_argument(
        "--bridge-times-s",
        type=float,
        nargs="+",
        default=list(DEFAULT_SOURCE_BRIDGE_TIMES_S),
    )
    bridge.add_argument("--out", type=str, default=None)
    bridge.set_defaults(func=command_casimir_source_bridge_report)

    return parser


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
