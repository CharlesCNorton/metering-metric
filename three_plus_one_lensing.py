"""
3+1D photon-sector lensing scaffold for static spherical metrics.

This file is the current theory-side workbench for the metering-metric
project's active lensing task. It evaluates null-geodesic observables for a
static spherically symmetric photon metric

    ds_gamma^2 = -A(r) dt^2 + B(r) dr^2 + r^2 dOmega^2

and makes one mainline point explicit:

- if photons decouple onto a flat bare metric, direct metering-induced photon
  deflection and gravitational redshift vanish
- if a nontrivial photon metric is proposed, its lensing and redshift can be
  computed here rather than inferred from lower-dimensional curvature alone
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import warnings
from dataclasses import dataclass
from pathlib import Path
from typing import Callable

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from scipy.stats import spearmanr
from scipy.integrate import IntegrationWarning, cumulative_trapezoid, quad
from scipy.ndimage import gaussian_filter, label, map_coordinates
from lensing import build_residual_map as archival_build_residual_map
from lensing import score_cutout as archival_score_cutout


MetricFn = Callable[[float], float]


GW170817_SPEED_LIMIT = 3.0e-15
VLBA_GAMMA_CENTRAL = 0.9998
VLBA_GAMMA_SIGMA = 3.0e-4
CASSINI_GAMMA_CENTRAL = 1.0 + 2.1e-5
CASSINI_GAMMA_SIGMA = 2.3e-5
GW170817_WEP_DELTA_GAMMA_MW = 5.9e-8
GW170817_WEP_DELTA_GAMMA_VIRGO = 0.9e-10
MICROSCOPE_ETA_CENTRAL = -1.5e-15
MICROSCOPE_ETA_SIGMA = math.sqrt((2.3e-15) ** 2 + (1.5e-15) ** 2)


@dataclass(frozen=True)
class StaticSphericalPhotonMetric:
    name: str
    lapse_A: MetricFn
    radial_B: MetricFn
    r_min: float
    asymptotic_lapse: float = 1.0
    metadata: dict[str, float] | None = None


@dataclass(frozen=True)
class SphericalMeteringProfile:
    radii: np.ndarray
    source: np.ndarray
    mu: np.ndarray
    mu_prime: np.ndarray
    screening_mass: float
    metadata: dict[str, float]


def tanh_lapse(mu: float, alpha: float, epsilon: float = 0.0) -> float:
    return float(epsilon + (1.0 - epsilon) * math.tanh(alpha * mu))


def gaussian_mu(r: float, mu0: float, sigma: float) -> float:
    return float(mu0 * math.exp(-(r * r) / (2.0 * sigma * sigma)))


def gaussian_source(r: np.ndarray, amplitude: float, sigma: float) -> np.ndarray:
    return amplitude * np.exp(-(r * r) / (2.0 * sigma * sigma))


def beta_model_source(r: np.ndarray, amplitude: float, core_radius: float, beta_index: float) -> np.ndarray:
    scale = np.maximum(core_radius, 1.0e-12)
    return amplitude * np.power(1.0 + (r / scale) ** 2, -1.5 * beta_index)


def softened_nfw_source(
    r: np.ndarray,
    amplitude: float,
    core_radius: float,
    scale_radius: float,
    outer_slope: float,
) -> np.ndarray:
    core = np.maximum(core_radius, 1.0e-12)
    scale = np.maximum(scale_radius, core)
    return amplitude / ((1.0 + (r / core)) * np.power(1.0 + (r / scale), outer_slope - 1.0))


def build_source_profile_values(
    radii: np.ndarray,
    source_kind: str,
    source_amplitude: float,
    source_sigma: float,
    source_core_radius: float,
    source_scale_radius: float,
    source_beta: float,
    source_outer_slope: float,
) -> np.ndarray:
    if source_kind == "gaussian":
        return gaussian_source(radii, amplitude=source_amplitude, sigma=source_sigma)
    if source_kind == "beta_model":
        return beta_model_source(radii, amplitude=source_amplitude, core_radius=source_core_radius, beta_index=source_beta)
    if source_kind == "softened_nfw":
        return softened_nfw_source(
            radii,
            amplitude=source_amplitude,
            core_radius=source_core_radius,
            scale_radius=source_scale_radius,
            outer_slope=source_outer_slope,
        )
    raise ValueError(f"Unknown source_kind: {source_kind}")


def reverse_cumulative_trapezoid(values: np.ndarray, grid: np.ndarray) -> np.ndarray:
    result = np.zeros_like(values)
    if values.size > 1:
        shell_steps = 0.5 * (values[1:] + values[:-1]) * np.diff(grid)
        result[:-1] = np.cumsum(shell_steps[::-1])[::-1]
    return result


def solve_screened_poisson_profile_from_source(
    source: np.ndarray,
    screening_mass: float,
    profile_r_max: float,
    profile_samples: int,
    source_metadata: dict[str, float | str],
) -> SphericalMeteringProfile:
    if profile_samples < 5:
        raise ValueError("profile_samples must be at least 5.")

    radii = np.linspace(0.0, profile_r_max, profile_samples, dtype=float)
    step = float(radii[1] - radii[0])
    if source.shape != radii.shape:
        raise ValueError("source must have the same shape as the radial grid.")

    interior_radii = radii[1:-1]
    interior_count = interior_radii.size
    diagonal = np.full(interior_count, (-2.0 / (step * step)) - (screening_mass * screening_mass), dtype=float)
    off_diagonal = np.full(interior_count - 1, 1.0 / (step * step), dtype=float)
    operator = np.diag(diagonal) + np.diag(off_diagonal, 1) + np.diag(off_diagonal, -1)
    rhs = -interior_radii * source[1:-1]

    u = np.zeros_like(radii)
    u[1:-1] = np.linalg.solve(operator, rhs)

    safe_radii = np.maximum(radii, 1.0e-12)
    mu = np.zeros_like(radii)
    mu[1:] = u[1:] / safe_radii[1:]
    mu[0] = mu[1]
    mu_prime = np.gradient(mu, step, edge_order=2)
    radial_laplacian = np.gradient(mu_prime, step, edge_order=2)
    radial_laplacian[1:] += 2.0 * mu_prime[1:] / safe_radii[1:]
    radial_laplacian[0] = radial_laplacian[1]
    residual = radial_laplacian - (screening_mass * screening_mass * mu) + source

    return SphericalMeteringProfile(
        radii=radii,
        source=source,
        mu=mu,
        mu_prime=mu_prime,
        screening_mass=screening_mass,
        metadata={
            "profile_r_max": float(profile_r_max),
            "profile_samples": float(profile_samples),
            "max_mu": float(np.max(mu)),
            "max_source": float(np.max(source)),
            "max_field_equation_residual": float(np.max(np.abs(residual))),
            **source_metadata,
        },
    )


def solve_screened_poisson_profile(
    source_amplitude: float,
    source_sigma: float,
    screening_mass: float,
    profile_r_max: float,
    profile_samples: int,
) -> SphericalMeteringProfile:
    radii = np.linspace(0.0, profile_r_max, profile_samples, dtype=float)
    source = gaussian_source(radii, amplitude=source_amplitude, sigma=source_sigma)
    return solve_screened_poisson_profile_from_source(
        source=source,
        screening_mass=screening_mass,
        profile_r_max=profile_r_max,
        profile_samples=profile_samples,
        source_metadata={
            "source_kind": "gaussian",
            "source_amplitude": float(source_amplitude),
            "source_sigma": float(source_sigma),
        },
    )


def solve_screened_poisson_cluster_profile(
    source_kind: str,
    source_amplitude: float,
    source_sigma: float,
    source_core_radius: float,
    source_scale_radius: float,
    source_beta: float,
    source_outer_slope: float,
    screening_mass: float,
    profile_r_max: float,
    profile_samples: int,
) -> SphericalMeteringProfile:
    radii = np.linspace(0.0, profile_r_max, profile_samples, dtype=float)
    source = build_source_profile_values(
        radii=radii,
        source_kind=source_kind,
        source_amplitude=source_amplitude,
        source_sigma=source_sigma,
        source_core_radius=source_core_radius,
        source_scale_radius=source_scale_radius,
        source_beta=source_beta,
        source_outer_slope=source_outer_slope,
    )
    metadata: dict[str, float | str] = {
        "source_kind": source_kind,
        "source_amplitude": float(source_amplitude),
        "source_sigma": float(source_sigma),
        "source_core_radius": float(source_core_radius),
        "source_scale_radius": float(source_scale_radius),
        "source_beta": float(source_beta),
        "source_outer_slope": float(source_outer_slope),
    }
    return solve_screened_poisson_profile_from_source(
        source=source,
        screening_mass=screening_mass,
        profile_r_max=profile_r_max,
        profile_samples=profile_samples,
        source_metadata=metadata,
    )


def flat_decoupled_metric() -> StaticSphericalPhotonMetric:
    return StaticSphericalPhotonMetric(
        name="flat_decoupled",
        lapse_A=lambda r: 1.0,
        radial_B=lambda r: 1.0,
        r_min=1.0e-6,
        asymptotic_lapse=1.0,
    )


def schwarzschild_metric(mass: float) -> StaticSphericalPhotonMetric:
    horizon = 2.0 * mass
    floor = horizon * 1.0001

    def lapse_A(r: float) -> float:
        return 1.0 - 2.0 * mass / r

    def radial_B(r: float) -> float:
        return 1.0 / lapse_A(r)

    return StaticSphericalPhotonMetric(
        name="schwarzschild",
        lapse_A=lapse_A,
        radial_B=radial_B,
        r_min=floor,
        asymptotic_lapse=1.0,
    )


def toy_photon_coupled_metering_metric(
    mu0: float,
    sigma: float,
    alpha: float,
    epsilon: float,
) -> StaticSphericalPhotonMetric:
    def lapse_A(r: float) -> float:
        lapse = tanh_lapse(gaussian_mu(r, mu0=mu0, sigma=sigma), alpha=alpha, epsilon=epsilon)
        return lapse * lapse

    def radial_B(r: float) -> float:
        lapse = tanh_lapse(gaussian_mu(r, mu0=mu0, sigma=sigma), alpha=alpha, epsilon=epsilon)
        return 1.0 / max(lapse * lapse, 1.0e-15)

    return StaticSphericalPhotonMetric(
        name="toy_photon_coupled_metering",
        lapse_A=lapse_A,
        radial_B=radial_B,
        r_min=1.0e-6,
        asymptotic_lapse=1.0,
    )


def gaussian_mu_prime(r: float, mu0: float, sigma: float) -> float:
    return float(-(r / (sigma * sigma)) * gaussian_mu(r, mu0=mu0, sigma=sigma))


def metering_energy_density(
    r: float,
    mu0: float,
    sigma: float,
    screening_mass: float,
    density_scale: float,
) -> float:
    mu_value = gaussian_mu(r, mu0=mu0, sigma=sigma)
    gradient = gaussian_mu_prime(r, mu0=mu0, sigma=sigma)
    return float(density_scale * (0.5 * gradient * gradient + 0.5 * screening_mass * screening_mass * mu_value * mu_value))


def spherical_newtonian_potential(
    radii: np.ndarray,
    density: np.ndarray,
    gravity_scale: float,
) -> tuple[np.ndarray, np.ndarray]:
    shell_density = 4.0 * math.pi * radii * radii * density
    enclosed_mass = cumulative_trapezoid(shell_density, radii, initial=0.0)
    tail_density = 4.0 * math.pi * radii * density
    tail_integral = np.zeros_like(radii)
    if radii.size > 1:
        shell_steps = 0.5 * (tail_density[1:] + tail_density[:-1]) * np.diff(radii)
        tail_integral[:-1] = np.cumsum(shell_steps[::-1])[::-1]
    safe_radii = np.maximum(radii, 1.0e-12)
    potential = -gravity_scale * (enclosed_mass / safe_radii + tail_integral)
    potential -= potential[-1]
    return potential, enclosed_mass


def toy_backreacted_metering_metric(
    mu0: float,
    sigma: float,
    screening_mass: float,
    density_scale: float,
    gravity_scale: float,
    profile_r_max: float,
    profile_samples: int,
) -> StaticSphericalPhotonMetric:
    r_min = 1.0e-4
    radii = np.geomspace(r_min, profile_r_max, profile_samples)
    density = np.asarray(
        [
            metering_energy_density(
                float(r),
                mu0=mu0,
                sigma=sigma,
                screening_mass=screening_mass,
                density_scale=density_scale,
            )
            for r in radii
        ],
        dtype=float,
    )
    potential, enclosed_mass = spherical_newtonian_potential(radii, density, gravity_scale=gravity_scale)
    lapse_values = np.maximum(1.0 + 2.0 * potential, 1.0e-8)
    radial_values = np.maximum(1.0 - 2.0 * potential, 1.0e-8)

    def lapse_A(r: float) -> float:
        return float(np.interp(r, radii, lapse_values, left=lapse_values[0], right=1.0))

    def radial_B(r: float) -> float:
        return float(np.interp(r, radii, radial_values, left=radial_values[0], right=1.0))

    total_mass = float(enclosed_mass[-1])
    minimum_lapse = float(np.min(np.sqrt(lapse_values)))
    scaled_total_mass = float(gravity_scale * total_mass)
    scaled_compactness = 2.0 * gravity_scale * enclosed_mass / np.maximum(radii, 1.0e-12)
    return StaticSphericalPhotonMetric(
        name="toy_backreacted_metering",
        lapse_A=lapse_A,
        radial_B=radial_B,
        r_min=r_min,
        asymptotic_lapse=1.0,
        metadata={
            "effective_mass_integral": total_mass,
            "gravity_scaled_mass": scaled_total_mass,
            "minimum_photon_lapse": minimum_lapse,
            "maximum_compactness_2GM_over_r": float(np.max(scaled_compactness)),
            "profile_r_max": float(profile_r_max),
        },
    )


def metering_stress_energy(
    profile: SphericalMeteringProfile,
    density_scale: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    potential = 0.5 * profile.screening_mass * profile.screening_mass * profile.mu * profile.mu
    kinetic = 0.5 * profile.mu_prime * profile.mu_prime
    rho = density_scale * (kinetic + potential)
    radial_pressure = density_scale * (kinetic - potential)
    tangential_pressure = density_scale * (-kinetic - potential)
    return rho, radial_pressure, tangential_pressure


def determine_exterior_matching_index(
    profile: SphericalMeteringProfile,
    rho: np.ndarray,
    mass_function: np.ndarray,
    tolerance: float,
) -> int:
    max_source = max(float(np.max(profile.source)), 1.0e-30)
    max_rho = max(float(np.max(rho)), 1.0e-30)
    total_mass = max(float(mass_function[-1]), 1.0e-30)
    source_ratio = profile.source / max_source
    rho_ratio = rho / max_rho
    tail_fraction = np.maximum(total_mass - mass_function, 0.0) / total_mass
    candidate = np.where(
        (source_ratio <= tolerance)
        & (rho_ratio <= tolerance)
        & (tail_fraction <= tolerance)
    )[0]
    if candidate.size == 0:
        return len(profile.radii) - 1
    first = int(candidate[0])
    return max(first, 2)


def einstein_scalar_backreaction_metric_from_profile(
    profile: SphericalMeteringProfile,
    density_scale: float,
    gravity_scale: float,
    exterior_match_tolerance: float = 1.0e-6,
) -> StaticSphericalPhotonMetric:
    radii = profile.radii
    safe_radii = np.maximum(radii, 1.0e-8)
    rho, radial_pressure, _ = metering_stress_energy(profile, density_scale=density_scale)
    source_density = 4.0 * math.pi * gravity_scale * radii * radii * rho
    mass_function = cumulative_trapezoid(source_density, radii, initial=0.0)
    compactness = 2.0 * mass_function / safe_radii
    compactness[0] = 0.0
    denominator = np.maximum(safe_radii * np.maximum(safe_radii - 2.0 * mass_function, 1.0e-10), 1.0e-10)
    phi_prime = (mass_function + 4.0 * math.pi * gravity_scale * radii * radii * radii * radial_pressure) / denominator
    phi_prime[0] = 0.0
    phi = -reverse_cumulative_trapezoid(phi_prime, radii)

    match_index = determine_exterior_matching_index(profile, rho=rho, mass_function=mass_function, tolerance=exterior_match_tolerance)
    match_radius = float(radii[match_index])
    match_mass = float(mass_function[match_index])
    exterior_lapse_match = max(1.0 - (2.0 * match_mass / max(match_radius, 1.0e-10)), 1.0e-12)
    phi = phi - phi[match_index] + 0.5 * math.log(exterior_lapse_match)

    lapse_values = np.maximum(np.exp(2.0 * phi), 1.0e-8)
    radial_values = 1.0 / np.maximum(1.0 - compactness, 1.0e-8)
    radial_values[0] = radial_values[1]

    interior_radii = radii[: match_index + 1]
    interior_lapse = lapse_values[: match_index + 1]
    interior_radial = radial_values[: match_index + 1]
    r_min = max(float(interior_radii[1]), 1.0e-6)

    def lapse_A(r: float) -> float:
        if r >= match_radius:
            return float(max(1.0 - (2.0 * match_mass / max(r, 1.0e-10)), 1.0e-12))
        return float(np.interp(r, interior_radii, interior_lapse, left=interior_lapse[1], right=interior_lapse[-1]))

    def radial_B(r: float) -> float:
        if r >= match_radius:
            denom = max(1.0 - (2.0 * match_mass / max(r, 1.0e-10)), 1.0e-12)
            return float(1.0 / denom)
        return float(np.interp(r, interior_radii, interior_radial, left=interior_radial[1], right=interior_radial[-1]))

    metadata = dict(profile.metadata)
    metadata.update(
        {
            "density_scale": float(density_scale),
            "gravity_scale": float(gravity_scale),
            "max_energy_density": float(np.max(rho)),
            "max_radial_pressure": float(np.max(radial_pressure)),
            "max_compactness_2GM_over_r": float(np.max(compactness)),
            "minimum_photon_lapse": float(np.min(np.sqrt(interior_lapse))),
            "gravity_scaled_mass": float(mass_function[-1]),
            "exterior_match_radius": match_radius,
            "exterior_match_mass": match_mass,
            "exterior_match_tolerance": float(exterior_match_tolerance),
        }
    )
    return StaticSphericalPhotonMetric(
        name="einstein_scalar_backreaction",
        lapse_A=lapse_A,
        radial_B=radial_B,
        r_min=r_min,
        asymptotic_lapse=1.0,
        metadata=metadata,
    )


def einstein_scalar_backreaction_metric(
    source_amplitude: float,
    source_sigma: float,
    screening_mass: float,
    density_scale: float,
    gravity_scale: float,
    profile_r_max: float,
    profile_samples: int,
) -> StaticSphericalPhotonMetric:
    profile = solve_screened_poisson_profile(
        source_amplitude=source_amplitude,
        source_sigma=source_sigma,
        screening_mass=screening_mass,
        profile_r_max=profile_r_max,
        profile_samples=profile_samples,
    )
    return einstein_scalar_backreaction_metric_from_profile(
        profile=profile,
        density_scale=density_scale,
        gravity_scale=gravity_scale,
    )


def photon_disformal_metering_metric(
    source_amplitude: float,
    source_sigma: float,
    screening_mass: float,
    temporal_coupling: float,
    radial_coupling: float,
    photon_alpha: float,
    profile_r_max: float,
    profile_samples: int,
) -> StaticSphericalPhotonMetric:
    profile = solve_screened_poisson_profile(
        source_amplitude=source_amplitude,
        source_sigma=source_sigma,
        screening_mass=screening_mass,
        profile_r_max=profile_r_max,
        profile_samples=profile_samples,
    )
    radii = profile.radii
    activation = np.tanh(photon_alpha * profile.mu)
    lapse_values = np.maximum(1.0 - temporal_coupling * activation * activation, 1.0e-8)
    radial_values = np.maximum(1.0 + radial_coupling * profile.mu_prime * profile.mu_prime, 1.0e-8)
    r_min = max(float(radii[1]), 1.0e-6)

    def lapse_A(r: float) -> float:
        return float(np.interp(r, radii, lapse_values, left=lapse_values[1], right=1.0))

    def radial_B(r: float) -> float:
        return float(np.interp(r, radii, radial_values, left=radial_values[1], right=1.0))

    metadata = dict(profile.metadata)
    metadata.update(
        {
            "temporal_coupling": float(temporal_coupling),
            "radial_coupling": float(radial_coupling),
            "photon_alpha": float(photon_alpha),
            "max_activation": float(np.max(activation * activation)),
            "minimum_photon_lapse": float(np.min(np.sqrt(lapse_values))),
            "maximum_photon_radial_factor": float(np.max(radial_values)),
        }
    )
    return StaticSphericalPhotonMetric(
        name="photon_disformal_metering",
        lapse_A=lapse_A,
        radial_B=radial_B,
        r_min=r_min,
        asymptotic_lapse=1.0,
        metadata=metadata,
    )


def radial_potential(metric: StaticSphericalPhotonMetric, r: float, impact_parameter: float) -> float:
    return (1.0 / metric.lapse_A(r)) - (impact_parameter * impact_parameter) / (r * r)


def find_turning_radius(
    metric: StaticSphericalPhotonMetric,
    impact_parameter: float,
    r_max: float,
    samples: int = 4096,
) -> float:
    lo = max(metric.r_min, 1.0e-6)
    grid = np.geomspace(lo, r_max, samples)
    values = np.asarray([radial_potential(metric, float(r), impact_parameter) for r in grid], dtype=float)
    finite = np.isfinite(values)
    grid = grid[finite]
    values = values[finite]
    if grid.size < 4:
        raise ValueError("Could not bracket a turning point.")

    sign_change_index = None
    for idx in range(len(values) - 1):
        left = float(values[idx])
        right = float(values[idx + 1])
        if left == 0.0:
            sign_change_index = idx
        elif left < 0.0 <= right:
            sign_change_index = idx
    if sign_change_index is None:
        raise ValueError("No photon turning point was found. The ray is captured or the bracket is too small.")

    a = float(grid[sign_change_index])
    b = float(grid[sign_change_index + 1])
    fa = radial_potential(metric, a, impact_parameter)
    fb = radial_potential(metric, b, impact_parameter)
    if fa == 0.0:
        return a
    if fb == 0.0:
        return b

    for _ in range(100):
        mid = 0.5 * (a + b)
        fm = radial_potential(metric, mid, impact_parameter)
        if not math.isfinite(fm):
            a = mid
            continue
        if abs(fm) < 1.0e-12:
            return mid
        if fa < 0.0 <= fm:
            b = mid
            fb = fm
        else:
            a = mid
            fa = fm
    return 0.5 * (a + b)


def deflection_integrand(metric: StaticSphericalPhotonMetric, r: float, impact_parameter: float) -> float:
    potential = radial_potential(metric, r, impact_parameter)
    if potential <= 0.0:
        return 0.0
    return impact_parameter * math.sqrt(metric.radial_B(r)) / (r * r * math.sqrt(potential))


def coordinate_time_integrand(metric: StaticSphericalPhotonMetric, r: float, impact_parameter: float) -> float:
    potential = radial_potential(metric, r, impact_parameter)
    if potential <= 0.0:
        return 0.0
    return math.sqrt(metric.radial_B(r)) / (metric.lapse_A(r) * math.sqrt(potential))


def integrate_from_turning_point(
    integrand: Callable[[float], float],
    turning_radius: float,
    upper_radius: float,
) -> tuple[float, float]:
    if upper_radius <= turning_radius:
        raise ValueError("upper_radius must exceed turning_radius.")
    upper_y = math.sqrt(upper_radius - turning_radius)

    def transformed(y: float) -> float:
        sample_y = max(y, 1.0e-12)
        radius = turning_radius + sample_y * sample_y
        return 2.0 * sample_y * integrand(radius)

    lower_y = min(1.0e-8, upper_y * 1.0e-8)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", IntegrationWarning)
        return quad(
            transformed,
            lower_y,
            upper_y,
            limit=400,
            epsabs=1.0e-10,
            epsrel=1.0e-10,
        )


def photon_deflection_angle(
    metric: StaticSphericalPhotonMetric,
    impact_parameter: float,
    r_max: float,
) -> dict:
    turning_radius = find_turning_radius(metric, impact_parameter=impact_parameter, r_max=r_max)
    metric_integral, metric_error = integrate_from_turning_point(
        lambda r: deflection_integrand(metric, r, impact_parameter),
        turning_radius=turning_radius,
        upper_radius=r_max,
    )
    flat_integral, flat_error = integrate_from_turning_point(
        lambda r: deflection_integrand(flat_decoupled_metric(), r, impact_parameter),
        turning_radius=impact_parameter,
        upper_radius=r_max,
    )
    angle = 2.0 * (metric_integral - flat_integral)
    return {
        "turning_radius": turning_radius,
        "deflection_angle_rad": angle,
        "deflection_angle_arcsec": angle * (180.0 / math.pi) * 3600.0,
        "integration_error_estimate": metric_error + flat_error,
    }


def photon_coordinate_time(
    metric: StaticSphericalPhotonMetric,
    impact_parameter: float,
    r_observer: float,
) -> dict:
    turning_radius = find_turning_radius(metric, impact_parameter=impact_parameter, r_max=r_observer)
    metric_time, metric_error = integrate_from_turning_point(
        lambda r: coordinate_time_integrand(metric, r, impact_parameter),
        turning_radius=turning_radius,
        upper_radius=r_observer,
    )
    flat_time, flat_error = integrate_from_turning_point(
        lambda r: coordinate_time_integrand(flat_decoupled_metric(), r, impact_parameter),
        turning_radius=impact_parameter,
        upper_radius=r_observer,
    )
    return {
        "turning_radius": turning_radius,
        "roundtrip_coordinate_time": 2.0 * metric_time,
        "flat_roundtrip_coordinate_time": 2.0 * flat_time,
        "excess_roundtrip_time": 2.0 * (metric_time - flat_time),
        "integration_error_estimate": metric_error + flat_error,
    }


def static_observer_redshift(metric: StaticSphericalPhotonMetric, r_emit: float, r_obs: float) -> dict:
    emit = metric.lapse_A(r_emit)
    obs = metric.lapse_A(r_obs)
    ratio = math.sqrt(emit / obs)
    return {
        "frequency_ratio_nu_obs_over_nu_emit": ratio,
        "redshift_z": (1.0 / ratio) - 1.0,
    }


def weak_field_schwarzschild_benchmark(mass: float, impact_parameter: float, r_max: float) -> dict:
    metric = schwarzschild_metric(mass)
    measured = photon_deflection_angle(metric, impact_parameter=impact_parameter, r_max=r_max)
    weak_field = 4.0 * mass / impact_parameter
    return {
        "metric": metric.name,
        "mass": mass,
        "impact_parameter": impact_parameter,
        "measured_deflection_angle_rad": measured["deflection_angle_rad"],
        "weak_field_prediction_rad": weak_field,
        "fractional_error_vs_4M_over_b": abs(measured["deflection_angle_rad"] - weak_field) / weak_field,
    }


def sample_axisymmetric_lensing_profile(
    metric: StaticSphericalPhotonMetric,
    impact_min: float,
    impact_max: float,
    samples: int,
    r_max: float,
) -> dict:
    if samples < 5:
        raise ValueError("samples must be at least 5.")
    impact_min = max(impact_min, metric.r_min * 1.01)
    impacts = np.geomspace(impact_min, impact_max, samples)
    deflections = np.zeros_like(impacts)
    delays = np.zeros_like(impacts)
    turning_radii = np.zeros_like(impacts)

    for idx, impact in enumerate(impacts):
        deflection = photon_deflection_angle(metric, impact_parameter=float(impact), r_max=r_max)
        timing = photon_coordinate_time(metric, impact_parameter=float(impact), r_observer=r_max)
        deflections[idx] = deflection["deflection_angle_rad"]
        delays[idx] = timing["excess_roundtrip_time"]
        turning_radii[idx] = deflection["turning_radius"]

    derivative = np.gradient(deflections, impacts, edge_order=2)
    convergence = 0.5 * (derivative + (deflections / impacts))
    shear_tangential = 0.5 * ((deflections / impacts) - derivative)
    magnification_inverse = (1.0 - convergence) ** 2 - shear_tangential**2
    threshold = max(1.0e-12, 1.0e-4 * float(np.max(np.abs(convergence))))
    convergence_resolved = np.where(np.abs(convergence) >= threshold, convergence, 0.0)
    sign_change_indices = np.where(np.signbit(convergence_resolved[:-1]) != np.signbit(convergence_resolved[1:]))[0]
    sign_flip_radii = [float(0.5 * (impacts[idx] + impacts[idx + 1])) for idx in sign_change_indices]

    return {
        "impact_parameters": impacts.tolist(),
        "turning_radii": turning_radii.tolist(),
        "deflection_angle_rad": deflections.tolist(),
        "excess_roundtrip_time": delays.tolist(),
        "convergence": convergence.tolist(),
        "shear_tangential": shear_tangential.tolist(),
        "magnification_inverse": magnification_inverse.tolist(),
        "sign_flip_radii": sign_flip_radii,
        "convergence_min": float(np.min(convergence)),
        "convergence_max": float(np.max(convergence)),
        "negative_convergence_fraction": float(np.mean(convergence < 0.0)),
        "resolved_zero_threshold": threshold,
        "negative_convergence_fraction_resolved": float(np.mean(convergence_resolved < 0.0)),
    }


def build_axisymmetric_lensing_map(
    radial_profile: dict,
    extent: float,
    grid_size: int,
    include_arrays: bool,
) -> dict:
    if grid_size < 5:
        raise ValueError("grid_size must be at least 5.")
    impacts = np.asarray(radial_profile["impact_parameters"], dtype=float)
    convergence = np.asarray(radial_profile["convergence"], dtype=float)
    shear_tangential = np.asarray(radial_profile["shear_tangential"], dtype=float)
    xs = np.linspace(-extent, extent, grid_size)
    ys = np.linspace(-extent, extent, grid_size)
    x_grid, y_grid = np.meshgrid(xs, ys)
    radius_grid = np.sqrt(x_grid * x_grid + y_grid * y_grid)
    kappa_grid = interpolate_radial_field(impacts, convergence, radius_grid)
    gamma_t_grid = interpolate_radial_field(impacts, shear_tangential, radius_grid)
    threshold = max(1.0e-12, 1.0e-4 * float(np.max(np.abs(convergence))))
    kappa_resolved = np.where(np.abs(kappa_grid) >= threshold, kappa_grid, 0.0)
    phi_grid = np.arctan2(y_grid, x_grid)
    gamma1_grid = -gamma_t_grid * np.cos(2.0 * phi_grid)
    gamma2_grid = -gamma_t_grid * np.sin(2.0 * phi_grid)

    result = {
        "extent": float(extent),
        "grid_size": int(grid_size),
        "kappa_min": float(np.min(kappa_grid)),
        "kappa_max": float(np.max(kappa_grid)),
        "gamma_t_max": float(np.max(np.abs(gamma_t_grid))),
        "negative_kappa_fraction": float(np.mean(kappa_grid < 0.0)),
        "resolved_zero_threshold": threshold,
        "negative_kappa_fraction_resolved": float(np.mean(kappa_resolved < 0.0)),
    }
    if include_arrays:
        result.update(
            {
                "x_coordinates": xs.tolist(),
                "y_coordinates": ys.tolist(),
                "kappa_map": kappa_grid.tolist(),
                "gamma1_map": gamma1_grid.tolist(),
                "gamma2_map": gamma2_grid.tolist(),
            }
        )
    return result


def build_radius_map(
    shape: tuple[int, int],
    center_x: float,
    center_y: float,
    pixel_scale: float,
) -> np.ndarray:
    yy, xx = np.indices(shape, dtype=float)
    x_coords = (xx - (shape[1] - 1) / 2.0) * pixel_scale
    y_coords = (yy - (shape[0] - 1) / 2.0) * pixel_scale
    return np.sqrt((x_coords - center_x) ** 2 + (y_coords - center_y) ** 2)


def gaussian_smooth_masked(image: np.ndarray, footprint: np.ndarray, sigma_px: float) -> np.ndarray:
    finite_mask = footprint & np.isfinite(image)
    if not np.any(finite_mask):
        raise ValueError("Masked Gaussian smoothing requires at least one finite pixel.")
    values = np.where(finite_mask, image, 0.0)
    weights = np.where(finite_mask, 1.0, 0.0)
    smooth_values = gaussian_filter(values, sigma=sigma_px)
    smooth_weights = gaussian_filter(weights, sigma=sigma_px)
    smooth = np.full(image.shape, np.nan, dtype=float)
    valid = smooth_weights > 1.0e-6
    smooth[valid] = smooth_values[valid] / smooth_weights[valid]
    smooth[~footprint] = np.nan
    return smooth


def build_residual_map(
    field: np.ndarray,
    footprint: np.ndarray,
    center_x: float,
    center_y: float,
    pixel_scale: float,
    radial_bin_width: float,
    smooth_sigma_px: float,
) -> tuple[np.ndarray, dict]:
    radius_map = build_radius_map(
        shape=field.shape,
        center_x=center_x,
        center_y=center_y,
        pixel_scale=pixel_scale,
    )
    finite_mask = footprint & np.isfinite(field) & np.isfinite(radius_map)
    if not np.any(finite_mask):
        raise ValueError("No finite pixels available for radial-median subtraction.")

    max_radius = float(np.max(radius_map[finite_mask]))
    bin_edges = np.arange(0.0, max_radius + radial_bin_width, radial_bin_width)
    if bin_edges.size < 2:
        bin_edges = np.array([0.0, max_radius + radial_bin_width], dtype=float)
    bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
    radial_medians = np.full(bin_centers.shape, np.nan, dtype=float)

    for index in range(bin_centers.size):
        in_bin = finite_mask & (radius_map >= bin_edges[index]) & (radius_map < bin_edges[index + 1])
        if np.any(in_bin):
            radial_medians[index] = float(np.median(field[in_bin]))

    valid = np.isfinite(radial_medians)
    if np.count_nonzero(valid) < 2:
        raise ValueError("Radial-median subtraction did not find enough populated bins.")

    baseline = np.interp(
        radius_map,
        bin_centers[valid],
        radial_medians[valid],
        left=float(radial_medians[valid][0]),
        right=float(radial_medians[valid][-1]),
    )
    radial_residual = field - baseline
    radial_residual[~footprint] = np.nan
    local_background = gaussian_smooth_masked(radial_residual, footprint, sigma_px=smooth_sigma_px)
    residual = radial_residual - local_background
    residual[~footprint] = np.nan
    return residual, {
        "radial_bin_width": float(radial_bin_width),
        "radial_bin_count": int(bin_centers.size),
        "populated_radial_bin_count": int(np.count_nonzero(valid)),
        "smooth_sigma_px": float(smooth_sigma_px),
    }


def interpolate_radial_field(impacts: np.ndarray, values: np.ndarray, radius: np.ndarray) -> np.ndarray:
    return np.interp(
        radius,
        impacts,
        values,
        left=float(values[0]),
        right=0.0,
    )


def component_weight(spec: dict) -> float:
    return float(spec.get("weight", spec.get("amplitude_scale", 1.0)))


def component_centroid(specs: list[dict]) -> tuple[float, float]:
    total_weight = max(sum(component_weight(spec) for spec in specs), 1.0e-12)
    center_x = sum(component_weight(spec) * float(spec["center_x"]) for spec in specs) / total_weight
    center_y = sum(component_weight(spec) * float(spec["center_y"]) for spec in specs) / total_weight
    return float(center_x), float(center_y)


def component_principal_axis_deg(specs: list[dict]) -> float:
    center_x, center_y = component_centroid(specs)
    total_weight = max(sum(component_weight(spec) for spec in specs), 1.0e-12)
    covariance = np.zeros((2, 2), dtype=float)
    for spec in specs:
        dx = float(spec["center_x"]) - center_x
        dy = float(spec["center_y"]) - center_y
        weight = component_weight(spec)
        displacement = np.array([[dx], [dy]], dtype=float)
        covariance += weight * (displacement @ displacement.T)
    covariance /= total_weight
    eigenvalues, eigenvectors = np.linalg.eigh(covariance)
    major_axis = eigenvectors[:, int(np.argmax(eigenvalues))]
    angle = math.degrees(math.atan2(float(major_axis[1]), float(major_axis[0])))
    return float(angle)


def residual_angular_harmonics(
    residual_map: np.ndarray,
    extent: float,
    center_x: float,
    center_y: float,
    max_order: int = 4,
) -> dict[str, float]:
    grid_size = residual_map.shape[0]
    pixel_scale = (2.0 * extent) / max(grid_size - 1, 1)
    yy, xx = np.indices(residual_map.shape, dtype=float)
    x_coords = (xx - (residual_map.shape[1] - 1) / 2.0) * pixel_scale - center_x
    y_coords = (yy - (residual_map.shape[0] - 1) / 2.0) * pixel_scale - center_y
    theta = np.arctan2(y_coords, x_coords)
    resolved_threshold = max(1.0e-12, 1.0e-4 * float(np.nanmax(np.abs(residual_map))))
    resolved = np.where(np.abs(residual_map) >= resolved_threshold, residual_map, 0.0)
    norm = float(np.sum(np.abs(resolved)))
    if norm <= 0.0:
        return {f"m{order}": 0.0 for order in range(1, max_order + 1)}

    harmonics: dict[str, float] = {}
    for order in range(1, max_order + 1):
        cosine = float(np.sum(resolved * np.cos(order * theta)) / norm)
        sine = float(np.sum(resolved * np.sin(order * theta)) / norm)
        harmonics[f"m{order}"] = float(math.sqrt(cosine * cosine + sine * sine))
    return harmonics


def count_connected_regions(mask: np.ndarray) -> int:
    structure = np.ones((3, 3), dtype=int)
    _, count = label(mask, structure=structure)
    return int(count)


def sample_field_along_axis(
    field: np.ndarray,
    extent: float,
    center_x: float,
    center_y: float,
    axis_angle_deg: float,
    radial_extent: float,
    samples: int,
) -> dict:
    if samples < 5:
        raise ValueError("samples must be at least 5.")
    offsets = np.linspace(-radial_extent, radial_extent, samples)
    angle = math.radians(axis_angle_deg)
    x_coords = center_x + offsets * math.cos(angle)
    y_coords = center_y + offsets * math.sin(angle)
    pixel_scale = (2.0 * extent) / max(field.shape[1] - 1, 1)
    x_index = (x_coords + extent) / pixel_scale
    y_index = (y_coords + extent) / pixel_scale
    samples_values = map_coordinates(field, [y_index, x_index], order=1, mode="nearest")
    threshold = max(1.0e-12, 1.0e-4 * float(np.nanmax(np.abs(field))))
    resolved = np.where(np.abs(samples_values) >= threshold, samples_values, 0.0)
    sign_changes = np.where(np.signbit(resolved[:-1]) != np.signbit(resolved[1:]))[0]
    sign_flip_positions = [float(0.5 * (offsets[idx] + offsets[idx + 1])) for idx in sign_changes]
    return {
        "axis_angle_deg": float(axis_angle_deg),
        "offsets": offsets.tolist(),
        "values": samples_values.tolist(),
        "resolved_zero_threshold": threshold,
        "sign_flip_positions": sign_flip_positions,
    }


def composite_component_specs(scenario: str, separation: float, subcluster_ratio: float) -> list[dict]:
    half = 0.5 * separation
    if scenario == "single_center":
        return [
            {"name": "core", "center_x": 0.0, "center_y": 0.0, "amplitude_scale": 1.0, "weight": 1.0},
        ]
    if scenario == "binary_equal":
        return [
            {"name": "west", "center_x": -half, "center_y": 0.0, "amplitude_scale": 1.0, "weight": 1.0},
            {"name": "east", "center_x": half, "center_y": 0.0, "amplitude_scale": 1.0, "weight": 1.0},
        ]
    if scenario == "bullet_offset":
        return [
            {"name": "main", "center_x": -0.25 * separation, "center_y": 0.0, "amplitude_scale": 1.0, "weight": 1.0},
            {
                "name": "subcluster",
                "center_x": 0.75 * separation,
                "center_y": 0.18 * separation,
                "amplitude_scale": subcluster_ratio,
                "sigma_scale": 0.75,
                "core_scale": 0.7,
                "scale_radius_scale": 0.7,
                "weight": subcluster_ratio,
            },
        ]
    if scenario == "triple_asymmetric":
        return [
            {"name": "core", "center_x": 0.0, "center_y": 0.0, "amplitude_scale": 1.0, "weight": 1.0},
            {
                "name": "northwest",
                "center_x": -0.7 * separation,
                "center_y": 0.45 * separation,
                "amplitude_scale": 0.55,
                "sigma_scale": 0.8,
                "core_scale": 0.8,
                "scale_radius_scale": 0.75,
                "weight": 0.55,
            },
            {
                "name": "southeast",
                "center_x": 0.9 * separation,
                "center_y": -0.35 * separation,
                "amplitude_scale": 0.4,
                "sigma_scale": 0.65,
                "core_scale": 0.6,
                "scale_radius_scale": 0.65,
                "weight": 0.4,
            },
        ]
    raise ValueError(f"Unknown scenario: {scenario}")


def load_member_geometry_component_specs(
    member_table_path: Path,
    top_members: int,
    mass_key: str,
    member_flag_key: str,
    ra_key: str,
    dec_key: str,
    amplitude_scaling_exponent: float,
    size_scaling_exponent: float,
    minimum_size_scale: float,
) -> tuple[list[dict], dict[str, float]]:
    with member_table_path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))

    members = []
    for row in rows:
        try:
            member_flag = int(float(row[member_flag_key]))
            mass = float(row[mass_key])
            ra = float(row[ra_key])
            dec = float(row[dec_key])
        except (KeyError, ValueError):
            continue
        if member_flag != 1:
            continue
        if not math.isfinite(mass) or mass <= 0.0:
            continue
        members.append(
            {
                "name": f"member_{row.get('IDcat', len(members) + 1)}",
                "mass": mass,
                "ra": ra,
                "dec": dec,
            }
        )

    if not members:
        raise ValueError("No valid cluster members were found in the member table.")

    members.sort(key=lambda item: item["mass"], reverse=True)
    selected = members[:top_members]
    total_mass = sum(item["mass"] for item in selected)
    centroid_ra = sum(item["mass"] * item["ra"] for item in selected) / total_mass
    centroid_dec = sum(item["mass"] * item["dec"] for item in selected) / total_mass
    cos_dec = math.cos(math.radians(centroid_dec))
    max_mass = max(item["mass"] for item in selected)

    component_specs = []
    for item in selected:
        weight = item["mass"] / max_mass
        amplitude_scale = weight**amplitude_scaling_exponent
        size_scale = max(weight**size_scaling_exponent, minimum_size_scale)
        component_specs.append(
            {
                "name": item["name"],
                "center_x": (item["ra"] - centroid_ra) * cos_dec * 3600.0,
                "center_y": (item["dec"] - centroid_dec) * 3600.0,
                "amplitude_scale": amplitude_scale,
                "sigma_scale": size_scale,
                "core_scale": size_scale,
                "scale_radius_scale": size_scale,
                "weight": weight,
                "member_mass": item["mass"],
                "layer": "member",
            }
        )

    metadata = {
        "top_members": float(len(component_specs)),
        "selected_mass_sum": float(total_mass),
        "max_member_mass": float(max_mass),
        "centroid_ra_deg": float(centroid_ra),
        "centroid_dec_deg": float(centroid_dec),
        "amplitude_scaling_exponent": float(amplitude_scaling_exponent),
        "size_scaling_exponent": float(size_scaling_exponent),
        "minimum_size_scale": float(minimum_size_scale),
    }
    return component_specs, metadata


def build_smoothed_background_component_specs(
    component_specs: list[dict],
    amplitude_fraction: float,
    size_multiplier: float,
    minimum_size_scale: float,
    position_shrink: float,
) -> list[dict]:
    if amplitude_fraction <= 0.0:
        return []

    clamped_position_shrink = min(max(position_shrink, 0.0), 1.0)
    background_specs = []
    for spec in component_specs:
        amplitude_scale = float(spec.get("amplitude_scale", 1.0)) * amplitude_fraction
        if amplitude_scale <= 0.0:
            continue
        size_scale = max(float(spec.get("sigma_scale", 1.0)) * size_multiplier, minimum_size_scale)
        background_specs.append(
            {
                **spec,
                "name": f"background_{spec['name']}",
                "center_x": (1.0 - clamped_position_shrink) * float(spec["center_x"]),
                "center_y": (1.0 - clamped_position_shrink) * float(spec["center_y"]),
                "amplitude_scale": amplitude_scale,
                "sigma_scale": size_scale,
                "core_scale": max(float(spec.get("core_scale", 1.0)) * size_multiplier, minimum_size_scale),
                "scale_radius_scale": max(
                    float(spec.get("scale_radius_scale", 1.0)) * size_multiplier,
                    minimum_size_scale,
                ),
                "weight": float(spec.get("weight", 1.0)) * amplitude_fraction,
                "layer": "smooth_background",
            }
        )
    return background_specs


def build_einstein_component_payloads(
    component_specs: list[dict],
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
) -> list[dict]:
    component_payloads = []
    for spec in component_specs:
        profile = solve_screened_poisson_cluster_profile(
            source_kind=source_kind,
            source_amplitude=source_amplitude * float(spec.get("amplitude_scale", 1.0)),
            source_sigma=source_sigma * float(spec.get("sigma_scale", 1.0)),
            source_core_radius=source_core_radius * float(spec.get("core_scale", 1.0)),
            source_scale_radius=source_scale_radius * float(spec.get("scale_radius_scale", 1.0)),
            source_beta=source_beta,
            source_outer_slope=source_outer_slope,
            screening_mass=screening_mass,
            profile_r_max=profile_r_max,
            profile_samples=profile_samples,
        )
        metric = einstein_scalar_backreaction_metric_from_profile(
            profile=profile,
            density_scale=density_scale,
            gravity_scale=gravity_scale,
            exterior_match_tolerance=exterior_match_tolerance,
        )
        radial_profile = sample_axisymmetric_lensing_profile(
            metric,
            impact_min=impact_min,
            impact_max=impact_max,
            samples=radial_samples,
            r_max=r_max,
        )
        component_payloads.append(
            {
                "spec": spec,
                "source_profile": profile.metadata,
                "metric_metadata": metric.metadata,
                "radial_profile": radial_profile,
            }
        )
    return component_payloads


def suggested_map_extent(component_specs: list[dict], padding_arcsec: float, minimum_extent: float) -> float:
    max_radius = max(
        math.hypot(float(spec["center_x"]), float(spec["center_y"]))
        for spec in component_specs
    )
    return float(max(minimum_extent, max_radius + padding_arcsec))


def build_multicomponent_lensing_map(
    component_payloads: list[dict],
    extent: float,
    grid_size: int,
    radial_bin_width: float,
    smooth_sigma_px: float,
    axis_samples: int,
    include_arrays: bool,
) -> dict:
    if grid_size < 5:
        raise ValueError("grid_size must be at least 5.")
    xs = np.linspace(-extent, extent, grid_size)
    ys = np.linspace(-extent, extent, grid_size)
    x_grid, y_grid = np.meshgrid(xs, ys)
    kappa_total = np.zeros_like(x_grid)
    gamma1_total = np.zeros_like(x_grid)
    gamma2_total = np.zeros_like(x_grid)

    for payload in component_payloads:
        radial_profile = payload["radial_profile"]
        impacts = np.asarray(radial_profile["impact_parameters"], dtype=float)
        convergence = np.asarray(radial_profile["convergence"], dtype=float)
        shear_tangential = np.asarray(radial_profile["shear_tangential"], dtype=float)
        dx = x_grid - float(payload["spec"]["center_x"])
        dy = y_grid - float(payload["spec"]["center_y"])
        radius = np.sqrt(dx * dx + dy * dy)
        kappa_local = interpolate_radial_field(impacts, convergence, radius)
        gamma_t_local = interpolate_radial_field(impacts, shear_tangential, radius)
        phi = np.arctan2(dy, dx)
        gamma1_local = -gamma_t_local * np.cos(2.0 * phi)
        gamma2_local = -gamma_t_local * np.sin(2.0 * phi)
        kappa_total += kappa_local
        gamma1_total += gamma1_local
        gamma2_total += gamma2_local

    gamma_t_total = np.sqrt(gamma1_total * gamma1_total + gamma2_total * gamma2_total)
    footprint = np.isfinite(kappa_total)
    pixel_scale = (2.0 * extent) / max(grid_size - 1, 1)
    center_x, center_y = component_centroid([payload["spec"] for payload in component_payloads])
    residual_map, residual_metadata = build_residual_map(
        field=kappa_total,
        footprint=footprint,
        center_x=center_x,
        center_y=center_y,
        pixel_scale=pixel_scale,
        radial_bin_width=radial_bin_width,
        smooth_sigma_px=smooth_sigma_px,
    )
    threshold = max(1.0e-12, 1.0e-4 * float(np.max(np.abs(kappa_total))))
    kappa_resolved = np.where(np.abs(kappa_total) >= threshold, kappa_total, 0.0)
    residual_threshold = max(1.0e-12, 1.0e-4 * float(np.nanmax(np.abs(residual_map))))
    residual_resolved = np.where(np.abs(residual_map) >= residual_threshold, residual_map, 0.0)
    axis_angle_deg = component_principal_axis_deg([payload["spec"] for payload in component_payloads])
    axis_profile = sample_field_along_axis(
        field=residual_map,
        extent=extent,
        center_x=center_x,
        center_y=center_y,
        axis_angle_deg=axis_angle_deg,
        radial_extent=extent,
        samples=axis_samples,
    )
    angular_harmonics = residual_angular_harmonics(
        residual_map=residual_map,
        extent=extent,
        center_x=center_x,
        center_y=center_y,
    )
    positive_regions = count_connected_regions(residual_resolved > 0.0)
    negative_regions = count_connected_regions(residual_resolved < 0.0)

    result = {
        "extent": float(extent),
        "grid_size": int(grid_size),
        "component_centroid": {"x": center_x, "y": center_y},
        "principal_axis_angle_deg": axis_angle_deg,
        "kappa_min": float(np.min(kappa_total)),
        "kappa_max": float(np.max(kappa_total)),
        "gamma_t_max": float(np.max(gamma_t_total)),
        "negative_kappa_fraction": float(np.mean(kappa_total < 0.0)),
        "negative_kappa_fraction_resolved": float(np.mean(kappa_resolved < 0.0)),
        "residual_min": float(np.nanmin(residual_map)),
        "residual_max": float(np.nanmax(residual_map)),
        "negative_residual_fraction": float(np.mean(residual_map < 0.0)),
        "negative_residual_fraction_resolved": float(np.mean(residual_resolved < 0.0)),
        "positive_residual_fraction_resolved": float(np.mean(residual_resolved > 0.0)),
        "residual_sign_flip_positions": axis_profile["sign_flip_positions"],
        "residual_axis_profile": axis_profile,
        "residual_angular_harmonics": angular_harmonics,
        "positive_resolved_region_count": positive_regions,
        "negative_resolved_region_count": negative_regions,
        "residual_map_metadata": residual_metadata,
        "resolved_zero_threshold": threshold,
        "resolved_residual_zero_threshold": residual_threshold,
    }
    if include_arrays:
        result.update(
            {
                "x_coordinates": xs.tolist(),
                "y_coordinates": ys.tolist(),
                "kappa_map": kappa_total.tolist(),
                "gamma1_map": gamma1_total.tolist(),
                "gamma2_map": gamma2_total.tolist(),
                "residual_map": residual_map.tolist(),
            }
        )
    return result


def sample_map_with_weighted_bilinear(
    image: np.ndarray,
    valid_mask: np.ndarray,
    x_coords: np.ndarray,
    y_coords: np.ndarray,
) -> np.ndarray:
    values = np.where(valid_mask, image, 0.0)
    weights = np.where(valid_mask, 1.0, 0.0)
    sampled_values = map_coordinates(values, [y_coords, x_coords], order=1, mode="constant", cval=0.0)
    sampled_weights = map_coordinates(weights, [y_coords, x_coords], order=1, mode="constant", cval=0.0)
    result = np.full(sampled_values.shape, np.nan, dtype=float)
    valid = sampled_weights > 1.0e-6
    result[valid] = sampled_values[valid] / sampled_weights[valid]
    return result


def rigid_transform_map(
    image: np.ndarray,
    extent: float,
    shift_x_arcsec: float,
    shift_y_arcsec: float,
    rotation_deg: float,
) -> np.ndarray:
    grid_size_y, grid_size_x = image.shape
    xs = np.linspace(-extent, extent, grid_size_x)
    ys = np.linspace(-extent, extent, grid_size_y)
    x_grid, y_grid = np.meshgrid(xs, ys)
    angle = math.radians(rotation_deg)
    cos_angle = math.cos(angle)
    sin_angle = math.sin(angle)
    x_shifted = x_grid - shift_x_arcsec
    y_shifted = y_grid - shift_y_arcsec
    source_x = cos_angle * x_shifted + sin_angle * y_shifted
    source_y = -sin_angle * x_shifted + cos_angle * y_shifted
    pixel_scale_x = (2.0 * extent) / max(grid_size_x - 1, 1)
    pixel_scale_y = (2.0 * extent) / max(grid_size_y - 1, 1)
    source_xpix = (source_x + extent) / pixel_scale_x
    source_ypix = (source_y + extent) / pixel_scale_y
    return sample_map_with_weighted_bilinear(
        image=image,
        valid_mask=np.isfinite(image),
        x_coords=source_xpix,
        y_coords=source_ypix,
    )


def build_archival_residual_map_on_member_grid(
    kappa_path: Path,
    centroid_ra_deg: float,
    centroid_dec_deg: float,
    extent: float,
    grid_size: int,
    residual_mode: str,
    smooth_sigma_px: float,
    radial_bin_arcsec: float,
) -> tuple[np.ndarray, dict[str, float], float]:
    with fits.open(kappa_path) as hdul:
        kappa = np.asarray(hdul[0].data, dtype=float)
        header = hdul[0].header
        wcs = WCS(header)

    pixel_scale_deg = header.get("CD1_1", header.get("CDELT1"))
    if pixel_scale_deg is None:
        raise ValueError("Could not determine pixel scale from FITS header.")
    pixel_scale_arcsec = abs(float(pixel_scale_deg)) * 3600.0
    footprint = np.isfinite(kappa)
    cluster_center_ra = float(header["CRVAL1"])
    cluster_center_dec = float(header["CRVAL2"])
    center_xpix, center_ypix = wcs.world_to_pixel_values(cluster_center_ra, cluster_center_dec)
    residual, residual_model = archival_build_residual_map(
        kappa=kappa,
        footprint=footprint,
        center_xpix=float(center_xpix),
        center_ypix=float(center_ypix),
        pixel_scale_arcsec=pixel_scale_arcsec,
        residual_mode=residual_mode,
        smooth_sigma_px=smooth_sigma_px,
        radial_bin_arcsec=radial_bin_arcsec,
    )

    xs = np.linspace(-extent, extent, grid_size)
    ys = np.linspace(-extent, extent, grid_size)
    x_grid, y_grid = np.meshgrid(xs, ys)
    cos_dec = math.cos(math.radians(centroid_dec_deg))
    ra_grid = centroid_ra_deg + x_grid / (3600.0 * max(cos_dec, 1.0e-8))
    dec_grid = centroid_dec_deg + y_grid / 3600.0
    xpix_grid, ypix_grid = wcs.world_to_pixel_values(ra_grid, dec_grid)
    archival_map = sample_map_with_weighted_bilinear(
        image=residual,
        valid_mask=np.isfinite(residual),
        x_coords=np.asarray(xpix_grid, dtype=float),
        y_coords=np.asarray(ypix_grid, dtype=float),
    )
    metadata: dict[str, float] = {
        "pixel_scale_arcsec": float(pixel_scale_arcsec),
        "residual_mean": float(np.nanmean(archival_map)),
        "residual_std": float(np.nanstd(archival_map)),
        "grid_extent_arcsec": float(extent),
        "grid_size": float(grid_size),
    }
    metadata.update({key: float(value) if isinstance(value, (int, float)) else value for key, value in residual_model.items()})
    return archival_map, metadata, pixel_scale_arcsec


def summarize_residual_map(
    residual_map: np.ndarray,
    extent: float,
    center_x: float,
    center_y: float,
) -> dict[str, float | list[float] | dict[str, float]]:
    residual_threshold = max(1.0e-12, 1.0e-4 * float(np.nanmax(np.abs(residual_map))))
    residual_resolved = np.where(np.abs(residual_map) >= residual_threshold, residual_map, 0.0)
    axis_profile = sample_field_along_axis(
        field=residual_map,
        extent=extent,
        center_x=center_x,
        center_y=center_y,
        axis_angle_deg=0.0,
        radial_extent=extent,
        samples=257,
    )
    return {
        "residual_min": float(np.nanmin(residual_map)),
        "residual_max": float(np.nanmax(residual_map)),
        "negative_residual_fraction_resolved": float(np.mean(residual_resolved < 0.0)),
        "positive_residual_fraction_resolved": float(np.mean(residual_resolved > 0.0)),
        "resolved_zero_threshold": float(residual_threshold),
        "residual_sign_flip_positions": axis_profile["sign_flip_positions"],
        "residual_angular_harmonics": residual_angular_harmonics(
            residual_map=residual_map,
            extent=extent,
            center_x=center_x,
            center_y=center_y,
        ),
        "positive_resolved_region_count": count_connected_regions(residual_resolved > 0.0),
        "negative_resolved_region_count": count_connected_regions(residual_resolved < 0.0),
    }


def safe_pearson_correlation(x: np.ndarray, y: np.ndarray) -> float:
    if x.size < 2 or y.size < 2:
        return float("nan")
    if float(np.std(x)) == 0.0 or float(np.std(y)) == 0.0:
        return float("nan")
    return float(np.corrcoef(x, y)[0, 1])


def safe_spearman_correlation(x: np.ndarray, y: np.ndarray) -> float:
    if x.size < 2 or y.size < 2:
        return float("nan")
    value = spearmanr(x, y, nan_policy="omit").statistic
    return float(value) if np.isfinite(value) else float("nan")


def linear_fit_summary(predictor: np.ndarray, response: np.ndarray) -> dict[str, float]:
    if predictor.size < 2 or response.size < 2:
        return {"intercept": float("nan"), "slope": float("nan"), "rmse": float("nan")}
    design = np.column_stack([np.ones(predictor.size), predictor])
    beta, *_ = np.linalg.lstsq(design, response, rcond=None)
    fitted = design @ beta
    rmse = float(np.sqrt(np.mean((response - fitted) ** 2)))
    return {"intercept": float(beta[0]), "slope": float(beta[1]), "rmse": rmse}


def jaccard_index(mask_a: np.ndarray, mask_b: np.ndarray) -> float:
    union = np.count_nonzero(mask_a | mask_b)
    if union == 0:
        return float("nan")
    return float(np.count_nonzero(mask_a & mask_b) / union)


def sample_member_scores_from_map(
    residual_map: np.ndarray,
    extent: float,
    member_specs: list[dict],
    inner_radius_arcsec: float,
    ring_inner_arcsec: float,
    ring_outer_arcsec: float,
) -> list[dict[str, float]]:
    pixel_scale = (2.0 * extent) / max(residual_map.shape[1] - 1, 1)
    footprint = np.isfinite(residual_map)
    rows = []
    for spec in member_specs:
        xpix = (float(spec["center_x"]) + extent) / pixel_scale
        ypix = (float(spec["center_y"]) + extent) / pixel_scale
        score = archival_score_cutout(
            residual=residual_map,
            footprint=footprint,
            xpix=xpix,
            ypix=ypix,
            pixel_scale_arcsec=pixel_scale,
            inner_radius_arcsec=inner_radius_arcsec,
            ring_inner_arcsec=ring_inner_arcsec,
            ring_outer_arcsec=ring_outer_arcsec,
        )
        if score is None:
            continue
        rows.append(
            {
                "name": str(spec["name"]),
                "member_mass": float(spec.get("member_mass", float("nan"))),
                "center_x": float(spec["center_x"]),
                "center_y": float(spec["center_y"]),
                **score,
            }
        )
    return rows


def compare_theory_and_archival_residuals(
    theory_map: np.ndarray,
    archival_map: np.ndarray,
    extent: float,
    member_specs: list[dict],
    inner_radius_arcsec: float,
    ring_inner_arcsec: float,
    ring_outer_arcsec: float,
) -> dict:
    finite_mask = np.isfinite(theory_map) & np.isfinite(archival_map)
    theory_values = theory_map[finite_mask]
    archival_values = archival_map[finite_mask]

    theory_threshold = max(1.0e-12, 1.0e-4 * float(np.nanmax(np.abs(theory_map))))
    archival_threshold = max(1.0e-12, 1.0e-4 * float(np.nanmax(np.abs(archival_map))))
    resolved_mask = finite_mask & (np.abs(theory_map) >= theory_threshold) & (np.abs(archival_map) >= archival_threshold)

    comparison = {
        "finite_pixel_count": int(np.count_nonzero(finite_mask)),
        "resolved_pixel_count": int(np.count_nonzero(resolved_mask)),
        "pearson_correlation": safe_pearson_correlation(theory_values, archival_values),
        "spearman_correlation": safe_spearman_correlation(theory_values, archival_values),
        "linear_fit_archive_from_theory": linear_fit_summary(theory_values, archival_values),
    }

    theory_norm = float(np.linalg.norm(theory_values))
    archival_norm = float(np.linalg.norm(archival_values))
    if theory_norm > 0.0 and archival_norm > 0.0:
        comparison["cosine_similarity"] = float(np.dot(theory_values, archival_values) / (theory_norm * archival_norm))
    else:
        comparison["cosine_similarity"] = float("nan")

    if np.count_nonzero(resolved_mask) > 0:
        theory_resolved = theory_map[resolved_mask]
        archival_resolved = archival_map[resolved_mask]
        comparison["resolved_sign_agreement_fraction"] = float(np.mean(np.signbit(theory_resolved) == np.signbit(archival_resolved)))
        comparison["positive_jaccard_resolved"] = jaccard_index(
            finite_mask & (theory_map >= theory_threshold),
            finite_mask & (archival_map >= archival_threshold),
        )
        comparison["negative_jaccard_resolved"] = jaccard_index(
            finite_mask & (theory_map <= -theory_threshold),
            finite_mask & (archival_map <= -archival_threshold),
        )
    else:
        comparison["resolved_sign_agreement_fraction"] = float("nan")
        comparison["positive_jaccard_resolved"] = float("nan")
        comparison["negative_jaccard_resolved"] = float("nan")

    theory_harmonics = residual_angular_harmonics(theory_map, extent=extent, center_x=0.0, center_y=0.0)
    archival_harmonics = residual_angular_harmonics(archival_map, extent=extent, center_x=0.0, center_y=0.0)
    theory_vector = np.asarray([theory_harmonics[f"m{order}"] for order in range(1, 5)], dtype=float)
    archival_vector = np.asarray([archival_harmonics[f"m{order}"] for order in range(1, 5)], dtype=float)
    comparison["harmonic_l2_distance"] = float(np.linalg.norm(theory_vector - archival_vector))
    comparison["theory_harmonics"] = theory_harmonics
    comparison["archival_harmonics"] = archival_harmonics

    theory_member_scores = sample_member_scores_from_map(
        residual_map=theory_map,
        extent=extent,
        member_specs=member_specs,
        inner_radius_arcsec=inner_radius_arcsec,
        ring_inner_arcsec=ring_inner_arcsec,
        ring_outer_arcsec=ring_outer_arcsec,
    )
    archival_member_scores = sample_member_scores_from_map(
        residual_map=archival_map,
        extent=extent,
        member_specs=member_specs,
        inner_radius_arcsec=inner_radius_arcsec,
        ring_inner_arcsec=ring_inner_arcsec,
        ring_outer_arcsec=ring_outer_arcsec,
    )
    archival_by_name = {row["name"]: row for row in archival_member_scores}
    paired = [
        {
            "name": row["name"],
            "member_mass": row["member_mass"],
            "theory_sign_flip_score": row["sign_flip_score"],
            "archival_sign_flip_score": archival_by_name[row["name"]]["sign_flip_score"],
            "theory_sign_flip_pass": row["sign_flip_pass"],
            "archival_sign_flip_pass": archival_by_name[row["name"]]["sign_flip_pass"],
        }
        for row in theory_member_scores
        if row["name"] in archival_by_name
    ]
    if paired:
        theory_scores = np.asarray([row["theory_sign_flip_score"] for row in paired], dtype=float)
        archival_scores = np.asarray([row["archival_sign_flip_score"] for row in paired], dtype=float)
        comparison["member_score_comparison"] = {
            "n_paired_members": int(len(paired)),
            "pearson_correlation": safe_pearson_correlation(theory_scores, archival_scores),
            "spearman_correlation": safe_spearman_correlation(theory_scores, archival_scores),
            "linear_fit_archive_from_theory": linear_fit_summary(theory_scores, archival_scores),
            "sign_flip_pass_agreement_fraction": float(
                np.mean(
                    np.asarray([row["theory_sign_flip_pass"] for row in paired], dtype=float)
                    == np.asarray([row["archival_sign_flip_pass"] for row in paired], dtype=float)
                )
            ),
            "paired_rows": paired,
        }
    else:
        comparison["member_score_comparison"] = {
            "n_paired_members": 0,
            "pearson_correlation": float("nan"),
            "spearman_correlation": float("nan"),
            "linear_fit_archive_from_theory": {"intercept": float("nan"), "slope": float("nan"), "rmse": float("nan")},
            "sign_flip_pass_agreement_fraction": float("nan"),
            "paired_rows": [],
        }

    return comparison


def comparison_objective(comparison: dict) -> float:
    positive_jaccard = float(comparison.get("positive_jaccard_resolved", float("nan")))
    negative_jaccard = float(comparison.get("negative_jaccard_resolved", float("nan")))
    if math.isfinite(positive_jaccard) and math.isfinite(negative_jaccard) and positive_jaccard >= 0.0 and negative_jaccard >= 0.0:
        overlap_geom = math.sqrt(positive_jaccard * negative_jaccard)
    else:
        overlap_geom = 0.0
    member_block = comparison.get("member_score_comparison", {})
    member_spearman = float(member_block.get("spearman_correlation", float("nan")))
    member_pass = float(member_block.get("sign_flip_pass_agreement_fraction", float("nan")))
    pixel_pearson = float(comparison.get("pearson_correlation", float("nan")))
    resolved_sign = float(comparison.get("resolved_sign_agreement_fraction", float("nan")))
    harmonic_l2 = float(comparison.get("harmonic_l2_distance", float("nan")))
    member_pass_term = member_pass if math.isfinite(member_pass) else 0.0
    member_spearman_term = 0.75 * max(member_spearman, 0.0) if math.isfinite(member_spearman) else 0.0
    resolved_sign_term = 0.40 * resolved_sign if math.isfinite(resolved_sign) else 0.0
    pixel_pearson_term = 0.15 * max(pixel_pearson, 0.0) if math.isfinite(pixel_pearson) else 0.0
    harmonic_penalty = 0.35 * harmonic_l2 if math.isfinite(harmonic_l2) else 0.0
    return float(
        member_pass_term
        + member_spearman_term
        + 0.50 * overlap_geom
        + resolved_sign_term
        + pixel_pearson_term
        - harmonic_penalty
    )


def optimize_rigid_alignment(
    theory_map: np.ndarray,
    archival_map: np.ndarray,
    extent: float,
    member_specs: list[dict],
    inner_radius_arcsec: float,
    ring_inner_arcsec: float,
    ring_outer_arcsec: float,
    max_shift_arcsec: float,
    shift_steps: int,
    max_rotation_deg: float,
    rotation_steps: int,
) -> dict:
    if max_shift_arcsec <= 0.0 and max_rotation_deg <= 0.0:
        return {
            "search_performed": False,
            "best_shift_x_arcsec": 0.0,
            "best_shift_y_arcsec": 0.0,
            "best_rotation_deg": 0.0,
            "best_objective": float("nan"),
            "best_comparison": {},
        }

    if shift_steps < 1:
        raise ValueError("shift_steps must be at least 1.")
    if rotation_steps < 1:
        raise ValueError("rotation_steps must be at least 1.")

    shift_values = np.linspace(-max_shift_arcsec, max_shift_arcsec, shift_steps)
    rotation_values = np.linspace(-max_rotation_deg, max_rotation_deg, rotation_steps)
    best_payload = None

    for rotation_deg in rotation_values:
        for shift_x_arcsec in shift_values:
            for shift_y_arcsec in shift_values:
                transformed = rigid_transform_map(
                    image=theory_map,
                    extent=extent,
                    shift_x_arcsec=float(shift_x_arcsec),
                    shift_y_arcsec=float(shift_y_arcsec),
                    rotation_deg=float(rotation_deg),
                )
                comparison = compare_theory_and_archival_residuals(
                    theory_map=transformed,
                    archival_map=archival_map,
                    extent=extent,
                    member_specs=member_specs,
                    inner_radius_arcsec=inner_radius_arcsec,
                    ring_inner_arcsec=ring_inner_arcsec,
                    ring_outer_arcsec=ring_outer_arcsec,
                )
                objective = comparison_objective(comparison)
                if best_payload is None or objective > best_payload["best_objective"]:
                    best_payload = {
                        "search_performed": True,
                        "best_shift_x_arcsec": float(shift_x_arcsec),
                        "best_shift_y_arcsec": float(shift_y_arcsec),
                        "best_rotation_deg": float(rotation_deg),
                        "best_objective": float(objective),
                        "best_comparison": comparison,
                    }

    assert best_payload is not None
    return best_payload


def asymptotic_speed_delta(metric: StaticSphericalPhotonMetric) -> dict:
    speed = math.sqrt(metric.asymptotic_lapse)
    delta = abs(speed - 1.0)
    return {
        "evaluation_radius": "asymptotic",
        "coordinate_speed": speed,
        "abs_speed_shift_from_c": delta,
        "gw170817_speed_limit": GW170817_SPEED_LIMIT,
        "passes_gw170817_speed": delta <= GW170817_SPEED_LIMIT,
    }


def effective_ppn_gamma(metric: StaticSphericalPhotonMetric, reference_radius: float) -> dict:
    lapse = metric.lapse_A(reference_radius)
    radial = metric.radial_B(reference_radius)
    delta_a = 1.0 - lapse
    delta_b = radial - 1.0
    if abs(delta_a) < 1.0e-18 and abs(delta_b) < 1.0e-18:
        gamma_eff = 1.0
        potential_u = 0.0
        mapped = True
    elif delta_a <= 0.0:
        gamma_eff = float("nan")
        potential_u = -0.5 * delta_a
        mapped = False
    else:
        gamma_eff = delta_b / delta_a
        potential_u = 0.5 * delta_a
        mapped = True

    payload = {
        "reference_radius": reference_radius,
        "lapse_A": lapse,
        "radial_B": radial,
        "effective_newtonian_potential_U": potential_u,
        "mapped_to_ppn_gamma": mapped,
        "gamma_eff": gamma_eff,
    }
    if mapped:
        delta_gamma = abs(gamma_eff - 1.0)
        payload.update(
            {
                "abs_delta_gamma_from_gr": delta_gamma,
                "vlba_gamma_central": VLBA_GAMMA_CENTRAL,
                "vlba_gamma_sigma": VLBA_GAMMA_SIGMA,
                "vlba_sigma_offset": abs(gamma_eff - VLBA_GAMMA_CENTRAL) / VLBA_GAMMA_SIGMA,
                "passes_vlba_2sigma": abs(gamma_eff - VLBA_GAMMA_CENTRAL) <= 2.0 * VLBA_GAMMA_SIGMA,
                "cassini_gamma_central": CASSINI_GAMMA_CENTRAL,
                "cassini_gamma_sigma": CASSINI_GAMMA_SIGMA,
                "cassini_sigma_offset": abs(gamma_eff - CASSINI_GAMMA_CENTRAL) / CASSINI_GAMMA_SIGMA,
                "passes_cassini_2sigma": abs(gamma_eff - CASSINI_GAMMA_CENTRAL) <= 2.0 * CASSINI_GAMMA_SIGMA,
                "gw170817_wep_mw_limit": GW170817_WEP_DELTA_GAMMA_MW,
                "passes_gw170817_wep_mw": delta_gamma <= GW170817_WEP_DELTA_GAMMA_MW,
                "gw170817_wep_virgo_limit": GW170817_WEP_DELTA_GAMMA_VIRGO,
                "passes_gw170817_wep_virgo": delta_gamma <= GW170817_WEP_DELTA_GAMMA_VIRGO,
            }
        )
    return payload


def photon_branch_required_radial_couplings(
    source_amplitude: float,
    source_sigma: float,
    screening_mass: float,
    temporal_coupling: float,
    photon_alpha: float,
    reference_radius: float,
    profile_r_max: float,
    profile_samples: int,
) -> dict:
    profile = solve_screened_poisson_profile(
        source_amplitude=source_amplitude,
        source_sigma=source_sigma,
        screening_mass=screening_mass,
        profile_r_max=profile_r_max,
        profile_samples=profile_samples,
    )
    mu_value = float(np.interp(reference_radius, profile.radii, profile.mu, left=profile.mu[1], right=profile.mu[-1]))
    mu_prime_value = float(np.interp(reference_radius, profile.radii, profile.mu_prime, left=profile.mu_prime[1], right=profile.mu_prime[-1]))
    activation_sq = math.tanh(photon_alpha * mu_value) ** 2
    slope_sq = mu_prime_value * mu_prime_value
    if slope_sq <= 0.0:
        return {
            "reference_radius": reference_radius,
            "mu_value": mu_value,
            "mu_prime": mu_prime_value,
            "activation_squared": activation_sq,
            "required_radial_coupling_gamma_gr": float("nan"),
        }

    def coupling_for_gamma(target_gamma: float) -> float:
        return target_gamma * temporal_coupling * activation_sq / slope_sq

    return {
        "reference_radius": reference_radius,
        "mu_value": mu_value,
        "mu_prime": mu_prime_value,
        "activation_squared": activation_sq,
        "required_radial_coupling_gamma_gr": coupling_for_gamma(1.0),
        "required_radial_coupling_vlba_central": coupling_for_gamma(VLBA_GAMMA_CENTRAL),
        "required_radial_coupling_vlba_minus_2sigma": coupling_for_gamma(VLBA_GAMMA_CENTRAL - 2.0 * VLBA_GAMMA_SIGMA),
        "required_radial_coupling_vlba_plus_2sigma": coupling_for_gamma(VLBA_GAMMA_CENTRAL + 2.0 * VLBA_GAMMA_SIGMA),
        "required_radial_coupling_cassini_central": coupling_for_gamma(CASSINI_GAMMA_CENTRAL),
        "required_radial_coupling_cassini_minus_2sigma": coupling_for_gamma(CASSINI_GAMMA_CENTRAL - 2.0 * CASSINI_GAMMA_SIGMA),
        "required_radial_coupling_cassini_plus_2sigma": coupling_for_gamma(CASSINI_GAMMA_CENTRAL + 2.0 * CASSINI_GAMMA_SIGMA),
        "required_radial_coupling_gw170817_wep_mw": coupling_for_gamma(1.0 + GW170817_WEP_DELTA_GAMMA_MW),
        "required_radial_coupling_gw170817_wep_virgo": coupling_for_gamma(1.0 + GW170817_WEP_DELTA_GAMMA_VIRGO),
    }


def photon_background_speed_bound(background_mu: float, temporal_coupling: float, photon_alpha: float) -> dict:
    activation_sq = math.tanh(photon_alpha * background_mu) ** 2
    lapse = max(1.0 - temporal_coupling * activation_sq, 1.0e-18)
    speed = math.sqrt(lapse)
    delta = abs(speed - 1.0)
    if activation_sq == 0.0:
        max_temporal = float("inf")
    else:
        max_temporal = (2.0 * GW170817_SPEED_LIMIT) / activation_sq
    return {
        "background_mu": background_mu,
        "activation_squared": activation_sq,
        "coordinate_speed": speed,
        "abs_speed_shift_from_c": delta,
        "gw170817_speed_limit": GW170817_SPEED_LIMIT,
        "passes_gw170817_speed": delta <= GW170817_SPEED_LIMIT,
        "max_temporal_coupling_from_gw170817": max_temporal,
    }


def scenario_metric(args: argparse.Namespace) -> StaticSphericalPhotonMetric:
    if args.scenario == "flat_decoupled":
        return flat_decoupled_metric()
    if args.scenario == "schwarzschild":
        return schwarzschild_metric(mass=args.mass)
    if args.scenario == "einstein_scalar_backreaction":
        return einstein_scalar_backreaction_metric(
            source_amplitude=args.source_amplitude,
            source_sigma=args.source_sigma,
            screening_mass=args.screening_mass,
            density_scale=args.density_scale,
            gravity_scale=args.gravity_scale,
            profile_r_max=args.profile_r_max,
            profile_samples=args.profile_samples,
        )
    if args.scenario == "photon_disformal_metering":
        return photon_disformal_metering_metric(
            source_amplitude=args.source_amplitude,
            source_sigma=args.source_sigma,
            screening_mass=args.screening_mass,
            temporal_coupling=args.temporal_coupling,
            radial_coupling=args.radial_coupling,
            photon_alpha=args.photon_alpha,
            profile_r_max=args.profile_r_max,
            profile_samples=args.profile_samples,
        )
    if args.scenario == "toy_backreacted_metering":
        return toy_backreacted_metering_metric(
            mu0=args.mu0,
            sigma=args.sigma,
            screening_mass=args.screening_mass,
            density_scale=args.density_scale,
            gravity_scale=args.gravity_scale,
            profile_r_max=args.profile_r_max,
            profile_samples=args.profile_samples,
        )
    if args.scenario == "toy_photon_coupled_metering":
        return toy_photon_coupled_metering_metric(
            mu0=args.mu0,
            sigma=args.sigma,
            alpha=args.alpha,
            epsilon=args.epsilon,
        )
    raise ValueError(f"Unknown scenario: {args.scenario}")


def command_evaluate(args: argparse.Namespace) -> None:
    metric = scenario_metric(args)
    result = {
        "scenario": metric.name,
        "impact_parameter": args.impact_parameter,
        "r_max": args.r_max,
        "r_emit": args.r_emit,
        "r_obs": args.r_obs,
        "deflection": photon_deflection_angle(metric, impact_parameter=args.impact_parameter, r_max=args.r_max),
        "coordinate_time": photon_coordinate_time(metric, impact_parameter=args.impact_parameter, r_observer=args.r_max),
        "redshift": static_observer_redshift(metric, r_emit=args.r_emit, r_obs=args.r_obs),
    }
    if metric.metadata is not None:
        result["metric_metadata"] = metric.metadata
    print(json.dumps(result, indent=2))


def command_benchmark(args: argparse.Namespace) -> None:
    result = weak_field_schwarzschild_benchmark(
        mass=args.mass,
        impact_parameter=args.impact_parameter,
        r_max=args.r_max,
    )
    print(json.dumps(result, indent=2))


def command_compare_completions(args: argparse.Namespace) -> None:
    flat_metric = flat_decoupled_metric()
    einstein_metric = einstein_scalar_backreaction_metric(
        source_amplitude=args.source_amplitude,
        source_sigma=args.source_sigma,
        screening_mass=args.screening_mass,
        density_scale=args.density_scale,
        gravity_scale=args.gravity_scale,
        profile_r_max=args.profile_r_max,
        profile_samples=args.profile_samples,
    )
    photon_metric = photon_disformal_metering_metric(
        source_amplitude=args.source_amplitude,
        source_sigma=args.source_sigma,
        screening_mass=args.screening_mass,
        temporal_coupling=args.temporal_coupling,
        radial_coupling=args.radial_coupling,
        photon_alpha=args.photon_alpha,
        profile_r_max=args.profile_r_max,
        profile_samples=args.profile_samples,
    )
    metrics = [flat_metric, einstein_metric, photon_metric]
    result = {
        "impact_parameter": args.impact_parameter,
        "r_max": args.r_max,
        "r_emit": args.r_emit,
        "r_obs": args.r_obs,
        "metrics": [],
    }
    for metric in metrics:
        payload = {
            "scenario": metric.name,
            "deflection": photon_deflection_angle(metric, impact_parameter=args.impact_parameter, r_max=args.r_max),
            "coordinate_time": photon_coordinate_time(metric, impact_parameter=args.impact_parameter, r_observer=args.r_max),
            "redshift": static_observer_redshift(metric, r_emit=args.r_emit, r_obs=args.r_obs),
        }
        if metric.metadata is not None:
            payload["metric_metadata"] = metric.metadata
        result["metrics"].append(payload)
    print(json.dumps(result, indent=2))


def command_constraint_audit(args: argparse.Namespace) -> None:
    reference_radius = args.reference_radius if args.reference_radius is not None else max(args.impact_parameter, args.source_sigma)
    flat_metric = flat_decoupled_metric()
    einstein_metric = einstein_scalar_backreaction_metric(
        source_amplitude=args.source_amplitude,
        source_sigma=args.source_sigma,
        screening_mass=args.screening_mass,
        density_scale=args.density_scale,
        gravity_scale=args.gravity_scale,
        profile_r_max=args.profile_r_max,
        profile_samples=args.profile_samples,
    )
    photon_metric = photon_disformal_metering_metric(
        source_amplitude=args.source_amplitude,
        source_sigma=args.source_sigma,
        screening_mass=args.screening_mass,
        temporal_coupling=args.temporal_coupling,
        radial_coupling=args.radial_coupling,
        photon_alpha=args.photon_alpha,
        profile_r_max=args.profile_r_max,
        profile_samples=args.profile_samples,
    )
    metrics = [flat_metric, einstein_metric, photon_metric]
    result = {
        "impact_parameter": args.impact_parameter,
        "reference_radius": reference_radius,
        "source_amplitude": args.source_amplitude,
            "source_sigma": args.source_sigma,
            "screening_mass": args.screening_mass,
            "density_scale": args.density_scale,
            "gravity_scale": args.gravity_scale,
            "temporal_coupling": args.temporal_coupling,
            "radial_coupling": args.radial_coupling,
            "photon_alpha": args.photon_alpha,
        "constraints": [],
        "photon_background_speed_bound": photon_background_speed_bound(
            background_mu=args.background_mu,
            temporal_coupling=args.temporal_coupling,
            photon_alpha=args.photon_alpha,
        ),
        "microscope_reference": {
            "eta_central": MICROSCOPE_ETA_CENTRAL,
            "eta_sigma_combined": MICROSCOPE_ETA_SIGMA,
            "note": "Current photon-side branches are not yet directly mapped to composition-dependent Eotvos ratios.",
        },
        "photon_branch_required_radial_couplings": photon_branch_required_radial_couplings(
            source_amplitude=args.source_amplitude,
            source_sigma=args.source_sigma,
            screening_mass=args.screening_mass,
            temporal_coupling=args.temporal_coupling,
            photon_alpha=args.photon_alpha,
            reference_radius=reference_radius,
            profile_r_max=args.profile_r_max,
            profile_samples=args.profile_samples,
        ),
    }
    for metric in metrics:
        payload = {
            "scenario": metric.name,
            "speed_bound": asymptotic_speed_delta(metric),
            "ppn_gamma_audit": effective_ppn_gamma(metric, reference_radius=reference_radius),
            "deflection": photon_deflection_angle(metric, impact_parameter=args.impact_parameter, r_max=args.r_max),
            "coordinate_time": photon_coordinate_time(metric, impact_parameter=args.impact_parameter, r_observer=args.r_max),
            "redshift": static_observer_redshift(metric, r_emit=args.r_emit, r_obs=args.r_obs),
        }
        if metric.metadata is not None:
            payload["metric_metadata"] = metric.metadata
        if metric.name == "einstein_scalar_backreaction":
            payload["interpretation_note"] = (
                "This branch is GR-derived. Outside the effective source tail it should approach a Schwarzschild exterior with gamma = 1. "
                "Local gamma offsets on the Gaussian shoulder are a profile/numerical extraction issue, not a distinct non-GR coupling."
            )
        if metric.name == "photon_disformal_metering":
            payload["interpretation_note"] = (
                "This branch is directly constrained by the gamma-based bounds. Surviving parameter space requires radial and temporal "
                "couplings to lie on a narrow tuning line set by the solved mu(r) profile at the chosen reference radius."
            )
        result["constraints"].append(payload)
    print(json.dumps(result, indent=2))


def command_cluster_map(args: argparse.Namespace) -> None:
    profile = solve_screened_poisson_cluster_profile(
        source_kind=args.source_kind,
        source_amplitude=args.source_amplitude,
        source_sigma=args.source_sigma,
        source_core_radius=args.source_core_radius,
        source_scale_radius=args.source_scale_radius,
        source_beta=args.source_beta,
        source_outer_slope=args.source_outer_slope,
        screening_mass=args.screening_mass,
        profile_r_max=args.profile_r_max,
        profile_samples=args.profile_samples,
    )
    metric = einstein_scalar_backreaction_metric_from_profile(
        profile=profile,
        density_scale=args.density_scale,
        gravity_scale=args.gravity_scale,
        exterior_match_tolerance=args.exterior_match_tolerance,
    )
    radial_profile = sample_axisymmetric_lensing_profile(
        metric,
        impact_min=args.impact_min,
        impact_max=args.impact_max,
        samples=args.radial_samples,
        r_max=args.r_max,
    )
    map_payload = build_axisymmetric_lensing_map(
        radial_profile=radial_profile,
        extent=args.map_extent,
        grid_size=args.grid_size,
        include_arrays=args.include_map_arrays,
    )
    result = {
        "scenario": "cluster_einstein_map",
        "source_profile": profile.metadata,
        "metric_metadata": metric.metadata,
        "radial_profile": radial_profile,
        "map": map_payload,
    }
    print(json.dumps(result, indent=2))


def command_cluster_composite_map(args: argparse.Namespace) -> None:
    component_specs = composite_component_specs(
        scenario=args.scenario,
        separation=args.separation,
        subcluster_ratio=args.subcluster_ratio,
    )
    component_payloads = []
    for spec in component_specs:
        profile = solve_screened_poisson_cluster_profile(
            source_kind=args.source_kind,
            source_amplitude=args.source_amplitude * float(spec.get("amplitude_scale", 1.0)),
            source_sigma=args.source_sigma * float(spec.get("sigma_scale", 1.0)),
            source_core_radius=args.source_core_radius * float(spec.get("core_scale", 1.0)),
            source_scale_radius=args.source_scale_radius * float(spec.get("scale_radius_scale", 1.0)),
            source_beta=args.source_beta,
            source_outer_slope=args.source_outer_slope,
            screening_mass=args.screening_mass,
            profile_r_max=args.profile_r_max,
            profile_samples=args.profile_samples,
        )
        metric = einstein_scalar_backreaction_metric_from_profile(
            profile=profile,
            density_scale=args.density_scale,
            gravity_scale=args.gravity_scale,
            exterior_match_tolerance=args.exterior_match_tolerance,
        )
        radial_profile = sample_axisymmetric_lensing_profile(
            metric,
            impact_min=args.impact_min,
            impact_max=args.impact_max,
            samples=args.radial_samples,
            r_max=args.r_max,
        )
        component_payloads.append(
            {
                "spec": spec,
                "source_profile": profile.metadata,
                "metric_metadata": metric.metadata,
                "radial_profile": radial_profile,
            }
        )

    map_payload = build_multicomponent_lensing_map(
        component_payloads=component_payloads,
        extent=args.map_extent,
        grid_size=args.grid_size,
        radial_bin_width=args.radial_bin_width,
        smooth_sigma_px=args.smooth_sigma_px,
        axis_samples=args.axis_samples,
        include_arrays=args.include_map_arrays,
    )
    result = {
        "scenario": f"composite_einstein_map:{args.scenario}",
        "source_kind": args.source_kind,
        "components": component_payloads,
        "map": map_payload,
    }
    print(json.dumps(result, indent=2))


def command_member_geometry_map(args: argparse.Namespace) -> None:
    component_specs, selection_metadata = load_member_geometry_component_specs(
        member_table_path=Path(args.member_table),
        top_members=args.top_members,
        mass_key=args.mass_key,
        member_flag_key=args.member_flag_key,
        ra_key=args.ra_key,
        dec_key=args.dec_key,
        amplitude_scaling_exponent=args.amplitude_scaling_exponent,
        size_scaling_exponent=args.size_scaling_exponent,
        minimum_size_scale=args.minimum_size_scale,
    )
    background_specs = build_smoothed_background_component_specs(
        component_specs=component_specs,
        amplitude_fraction=args.smooth_background_amplitude_fraction,
        size_multiplier=args.smooth_background_size_multiplier,
        minimum_size_scale=args.smooth_background_minimum_size_scale,
        position_shrink=args.smooth_background_position_shrink,
    )
    all_component_specs = [*component_specs, *background_specs]
    component_payloads = build_einstein_component_payloads(
        component_specs=all_component_specs,
        source_kind=args.source_kind,
        source_amplitude=args.source_amplitude,
        source_sigma=args.source_sigma,
        source_core_radius=args.source_core_radius,
        source_scale_radius=args.source_scale_radius,
        source_beta=args.source_beta,
        source_outer_slope=args.source_outer_slope,
        screening_mass=args.screening_mass,
        density_scale=args.density_scale,
        gravity_scale=args.gravity_scale,
        exterior_match_tolerance=args.exterior_match_tolerance,
        profile_r_max=args.profile_r_max,
        profile_samples=args.profile_samples,
        impact_min=args.impact_min,
        impact_max=args.impact_max,
        radial_samples=args.radial_samples,
        r_max=args.r_max,
    )

    map_payload = build_multicomponent_lensing_map(
        component_payloads=component_payloads,
        extent=args.map_extent,
        grid_size=args.grid_size,
        radial_bin_width=args.radial_bin_width,
        smooth_sigma_px=args.smooth_sigma_px,
        axis_samples=args.axis_samples,
        include_arrays=args.include_map_arrays,
    )
    result = {
        "scenario": "hff_member_geometry",
        "member_table": str(args.member_table),
        "selection_metadata": selection_metadata,
        "smooth_background": {
            "amplitude_fraction": args.smooth_background_amplitude_fraction,
            "size_multiplier": args.smooth_background_size_multiplier,
            "minimum_size_scale": args.smooth_background_minimum_size_scale,
            "position_shrink": args.smooth_background_position_shrink,
            "n_components": len(background_specs),
        },
        "components": component_payloads,
        "map": map_payload,
    }
    print(json.dumps(result, indent=2))


def command_compare_hff_member_geometry(args: argparse.Namespace) -> None:
    component_specs, selection_metadata = load_member_geometry_component_specs(
        member_table_path=Path(args.member_table),
        top_members=args.top_members,
        mass_key=args.mass_key,
        member_flag_key=args.member_flag_key,
        ra_key=args.ra_key,
        dec_key=args.dec_key,
        amplitude_scaling_exponent=args.amplitude_scaling_exponent,
        size_scaling_exponent=args.size_scaling_exponent,
        minimum_size_scale=args.minimum_size_scale,
    )
    background_specs = build_smoothed_background_component_specs(
        component_specs=component_specs,
        amplitude_fraction=args.smooth_background_amplitude_fraction,
        size_multiplier=args.smooth_background_size_multiplier,
        minimum_size_scale=args.smooth_background_minimum_size_scale,
        position_shrink=args.smooth_background_position_shrink,
    )
    all_component_specs = [*component_specs, *background_specs]
    map_extent = float(args.map_extent)
    if map_extent <= 0.0:
        map_extent = suggested_map_extent(
            component_specs=all_component_specs,
            padding_arcsec=args.auto_extent_padding_arcsec,
            minimum_extent=args.minimum_extent_arcsec,
        )
    component_payloads = build_einstein_component_payloads(
        component_specs=all_component_specs,
        source_kind=args.source_kind,
        source_amplitude=args.source_amplitude,
        source_sigma=args.source_sigma,
        source_core_radius=args.source_core_radius,
        source_scale_radius=args.source_scale_radius,
        source_beta=args.source_beta,
        source_outer_slope=args.source_outer_slope,
        screening_mass=args.screening_mass,
        density_scale=args.density_scale,
        gravity_scale=args.gravity_scale,
        exterior_match_tolerance=args.exterior_match_tolerance,
        profile_r_max=args.profile_r_max,
        profile_samples=args.profile_samples,
        impact_min=args.impact_min,
        impact_max=args.impact_max,
        radial_samples=args.radial_samples,
        r_max=args.r_max,
    )

    theory_map_payload = build_multicomponent_lensing_map(
        component_payloads=component_payloads,
        extent=map_extent,
        grid_size=args.grid_size,
        radial_bin_width=args.radial_bin_width,
        smooth_sigma_px=args.smooth_sigma_px,
        axis_samples=args.axis_samples,
        include_arrays=args.include_map_arrays,
    )
    theory_map = np.asarray(theory_map_payload["residual_map"] if args.include_map_arrays else build_multicomponent_lensing_map(
        component_payloads=component_payloads,
        extent=map_extent,
        grid_size=args.grid_size,
        radial_bin_width=args.radial_bin_width,
        smooth_sigma_px=args.smooth_sigma_px,
        axis_samples=args.axis_samples,
        include_arrays=True,
    )["residual_map"], dtype=float)

    archival_map, archival_metadata, archival_pixel_scale = build_archival_residual_map_on_member_grid(
        kappa_path=Path(args.kappa_map),
        centroid_ra_deg=float(selection_metadata["centroid_ra_deg"]),
        centroid_dec_deg=float(selection_metadata["centroid_dec_deg"]),
        extent=map_extent,
        grid_size=args.grid_size,
        residual_mode=args.residual_mode,
        smooth_sigma_px=args.archive_smooth_sigma_px,
        radial_bin_arcsec=args.archive_radial_bin_arcsec,
    )
    archival_summary = summarize_residual_map(
        residual_map=archival_map,
        extent=map_extent,
        center_x=0.0,
        center_y=0.0,
    )
    comparison = compare_theory_and_archival_residuals(
        theory_map=theory_map,
        archival_map=archival_map,
        extent=map_extent,
        member_specs=component_specs,
        inner_radius_arcsec=args.inner_radius_arcsec,
        ring_inner_arcsec=args.ring_inner_arcsec,
        ring_outer_arcsec=args.ring_outer_arcsec,
    )
    alignment = optimize_rigid_alignment(
        theory_map=theory_map,
        archival_map=archival_map,
        extent=map_extent,
        member_specs=component_specs,
        inner_radius_arcsec=args.inner_radius_arcsec,
        ring_inner_arcsec=args.ring_inner_arcsec,
        ring_outer_arcsec=args.ring_outer_arcsec,
        max_shift_arcsec=args.alignment_max_shift_arcsec,
        shift_steps=args.alignment_shift_steps,
        max_rotation_deg=args.alignment_max_rotation_deg,
        rotation_steps=args.alignment_rotation_steps,
    )

    result = {
        "scenario": "compare_hff_member_geometry",
        "cluster_label": args.cluster_label,
        "member_table": str(args.member_table),
        "kappa_map": str(args.kappa_map),
        "selection_metadata": selection_metadata,
        "parameters": {
            "top_members": args.top_members,
            "smooth_background_amplitude_fraction": args.smooth_background_amplitude_fraction,
            "smooth_background_size_multiplier": args.smooth_background_size_multiplier,
            "smooth_background_minimum_size_scale": args.smooth_background_minimum_size_scale,
            "smooth_background_position_shrink": args.smooth_background_position_shrink,
            "smooth_background_component_count": len(background_specs),
            "map_extent_arcsec": map_extent,
            "grid_size": args.grid_size,
            "source_kind": args.source_kind,
            "archive_residual_mode": args.residual_mode,
            "archive_smooth_sigma_px": args.archive_smooth_sigma_px,
            "archive_radial_bin_arcsec": args.archive_radial_bin_arcsec,
            "score_inner_radius_arcsec": args.inner_radius_arcsec,
            "score_ring_inner_arcsec": args.ring_inner_arcsec,
            "score_ring_outer_arcsec": args.ring_outer_arcsec,
            "alignment_max_shift_arcsec": args.alignment_max_shift_arcsec,
            "alignment_shift_steps": args.alignment_shift_steps,
            "alignment_max_rotation_deg": args.alignment_max_rotation_deg,
            "alignment_rotation_steps": args.alignment_rotation_steps,
        },
        "theory_map": theory_map_payload,
        "archival_map": {
            "metadata": archival_metadata,
            "pixel_scale_arcsec": archival_pixel_scale,
            **archival_summary,
        },
        "comparison": comparison,
        "alignment_search": alignment,
    }
    if args.include_map_arrays:
        result["archival_map"]["residual_map"] = archival_map.tolist()
    print(json.dumps(result, indent=2))


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Evaluate 3+1D photon observables for static spherical metrics.",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    evaluate = subparsers.add_parser("evaluate", help="Evaluate deflection, time delay, and redshift for a metric.")
    evaluate.add_argument(
        "--scenario",
        choices=[
            "flat_decoupled",
            "schwarzschild",
            "toy_photon_coupled_metering",
            "toy_backreacted_metering",
            "einstein_scalar_backreaction",
            "photon_disformal_metering",
        ],
        required=True,
    )
    evaluate.add_argument("--impact-parameter", type=float, default=50.0)
    evaluate.add_argument("--r-max", type=float, default=5.0e4)
    evaluate.add_argument("--r-emit", type=float, default=10.0)
    evaluate.add_argument("--r-obs", type=float, default=1.0e6)
    evaluate.add_argument("--mass", type=float, default=1.0)
    evaluate.add_argument("--mu0", type=float, default=1.0)
    evaluate.add_argument("--sigma", type=float, default=100.0)
    evaluate.add_argument("--alpha", type=float, default=5.0)
    evaluate.add_argument("--epsilon", type=float, default=0.05)
    evaluate.add_argument("--source-amplitude", type=float, default=1.0e-3)
    evaluate.add_argument("--source-sigma", type=float, default=120.0)
    evaluate.add_argument("--screening-mass", type=float, default=0.05)
    evaluate.add_argument("--density-scale", type=float, default=1.0e-6)
    evaluate.add_argument("--gravity-scale", type=float, default=1.0)
    evaluate.add_argument("--temporal-coupling", type=float, default=0.02)
    evaluate.add_argument("--radial-coupling", type=float, default=0.1)
    evaluate.add_argument("--photon-alpha", type=float, default=1.0)
    evaluate.add_argument("--profile-r-max", type=float, default=2.0e3)
    evaluate.add_argument("--profile-samples", type=int, default=4096)
    evaluate.set_defaults(func=command_evaluate)

    benchmark = subparsers.add_parser("benchmark-schwarzschild", help="Check the null-geodesic integral against 4M/b.")
    benchmark.add_argument("--mass", type=float, default=1.0)
    benchmark.add_argument("--impact-parameter", type=float, default=1000.0)
    benchmark.add_argument("--r-max", type=float, default=1.0e6)
    benchmark.set_defaults(func=command_benchmark)

    compare = subparsers.add_parser("compare-completions", help="Compare flat, Einstein-backreacted, and explicit photon-coupled branches.")
    compare.add_argument("--impact-parameter", type=float, default=50.0)
    compare.add_argument("--r-max", type=float, default=5.0e4)
    compare.add_argument("--r-emit", type=float, default=10.0)
    compare.add_argument("--r-obs", type=float, default=1.0e6)
    compare.add_argument("--source-amplitude", type=float, default=1.0e-3)
    compare.add_argument("--source-sigma", type=float, default=120.0)
    compare.add_argument("--screening-mass", type=float, default=0.05)
    compare.add_argument("--density-scale", type=float, default=1.0e-6)
    compare.add_argument("--gravity-scale", type=float, default=1.0)
    compare.add_argument("--temporal-coupling", type=float, default=0.02)
    compare.add_argument("--radial-coupling", type=float, default=0.1)
    compare.add_argument("--photon-alpha", type=float, default=1.0)
    compare.add_argument("--profile-r-max", type=float, default=2.0e3)
    compare.add_argument("--profile-samples", type=int, default=4096)
    compare.set_defaults(func=command_compare_completions)

    audit = subparsers.add_parser("constraint-audit", help="Audit the branches against GW170817, VLBA, Cassini, and WEP-style bounds.")
    audit.add_argument("--impact-parameter", type=float, default=120.0)
    audit.add_argument("--r-max", type=float, default=5.0e4)
    audit.add_argument("--r-emit", type=float, default=120.0)
    audit.add_argument("--r-obs", type=float, default=5.0e4)
    audit.add_argument("--reference-radius", type=float, default=None)
    audit.add_argument("--source-amplitude", type=float, default=1.0e-3)
    audit.add_argument("--source-sigma", type=float, default=120.0)
    audit.add_argument("--screening-mass", type=float, default=0.05)
    audit.add_argument("--density-scale", type=float, default=1.0e-6)
    audit.add_argument("--gravity-scale", type=float, default=1.0)
    audit.add_argument("--temporal-coupling", type=float, default=0.02)
    audit.add_argument("--radial-coupling", type=float, default=0.1)
    audit.add_argument("--photon-alpha", type=float, default=1.0)
    audit.add_argument("--background-mu", type=float, default=0.0)
    audit.add_argument("--profile-r-max", type=float, default=5.0e3)
    audit.add_argument("--profile-samples", type=int, default=4096)
    audit.set_defaults(func=command_constraint_audit)

    cluster_map = subparsers.add_parser("cluster-map", help="Build axisymmetric Einstein-branch lensing maps for cluster-style source profiles.")
    cluster_map.add_argument("--source-kind", choices=["gaussian", "beta_model", "softened_nfw"], default="softened_nfw")
    cluster_map.add_argument("--source-amplitude", type=float, default=1.0e-3)
    cluster_map.add_argument("--source-sigma", type=float, default=120.0)
    cluster_map.add_argument("--source-core-radius", type=float, default=60.0)
    cluster_map.add_argument("--source-scale-radius", type=float, default=350.0)
    cluster_map.add_argument("--source-beta", type=float, default=0.8)
    cluster_map.add_argument("--source-outer-slope", type=float, default=3.0)
    cluster_map.add_argument("--screening-mass", type=float, default=0.01)
    cluster_map.add_argument("--density-scale", type=float, default=1.0e-6)
    cluster_map.add_argument("--gravity-scale", type=float, default=1.0)
    cluster_map.add_argument("--exterior-match-tolerance", type=float, default=1.0e-6)
    cluster_map.add_argument("--profile-r-max", type=float, default=8.0e3)
    cluster_map.add_argument("--profile-samples", type=int, default=4096)
    cluster_map.add_argument("--impact-min", type=float, default=30.0)
    cluster_map.add_argument("--impact-max", type=float, default=4.0e3)
    cluster_map.add_argument("--radial-samples", type=int, default=64)
    cluster_map.add_argument("--r-max", type=float, default=8.0e4)
    cluster_map.add_argument("--map-extent", type=float, default=4.0e3)
    cluster_map.add_argument("--grid-size", type=int, default=65)
    cluster_map.add_argument("--include-map-arrays", action="store_true")
    cluster_map.set_defaults(func=command_cluster_map)

    composite_map = subparsers.add_parser(
        "cluster-composite-map",
        help="Build multi-component Einstein-branch residual maps for non-spherical cluster geometries.",
    )
    composite_map.add_argument(
        "--scenario",
        choices=["single_center", "binary_equal", "bullet_offset", "triple_asymmetric"],
        default="bullet_offset",
    )
    composite_map.add_argument("--separation", type=float, default=850.0)
    composite_map.add_argument("--subcluster-ratio", type=float, default=0.45)
    composite_map.add_argument("--source-kind", choices=["gaussian", "beta_model", "softened_nfw"], default="softened_nfw")
    composite_map.add_argument("--source-amplitude", type=float, default=1.0e-3)
    composite_map.add_argument("--source-sigma", type=float, default=120.0)
    composite_map.add_argument("--source-core-radius", type=float, default=60.0)
    composite_map.add_argument("--source-scale-radius", type=float, default=350.0)
    composite_map.add_argument("--source-beta", type=float, default=0.8)
    composite_map.add_argument("--source-outer-slope", type=float, default=3.0)
    composite_map.add_argument("--screening-mass", type=float, default=0.01)
    composite_map.add_argument("--density-scale", type=float, default=1.0e-6)
    composite_map.add_argument("--gravity-scale", type=float, default=1.0)
    composite_map.add_argument("--exterior-match-tolerance", type=float, default=1.0e-6)
    composite_map.add_argument("--profile-r-max", type=float, default=8.0e3)
    composite_map.add_argument("--profile-samples", type=int, default=4096)
    composite_map.add_argument("--impact-min", type=float, default=30.0)
    composite_map.add_argument("--impact-max", type=float, default=4.0e3)
    composite_map.add_argument("--radial-samples", type=int, default=64)
    composite_map.add_argument("--r-max", type=float, default=8.0e4)
    composite_map.add_argument("--map-extent", type=float, default=4.0e3)
    composite_map.add_argument("--grid-size", type=int, default=129)
    composite_map.add_argument("--radial-bin-width", type=float, default=120.0)
    composite_map.add_argument("--smooth-sigma-px", type=float, default=3.0)
    composite_map.add_argument("--axis-samples", type=int, default=257)
    composite_map.add_argument("--include-map-arrays", action="store_true")
    composite_map.set_defaults(func=command_cluster_composite_map)

    member_geometry = subparsers.add_parser(
        "member-geometry-map",
        help="Build Einstein-branch residual maps from a real cluster-member table.",
    )
    member_geometry.add_argument("--member-table", required=True)
    member_geometry.add_argument("--top-members", type=int, default=25)
    member_geometry.add_argument("--mass-key", default="mass_neb")
    member_geometry.add_argument("--member-flag-key", default="is_cluster_member")
    member_geometry.add_argument("--ra-key", default="ra")
    member_geometry.add_argument("--dec-key", default="dec")
    member_geometry.add_argument("--amplitude-scaling-exponent", type=float, default=1.0)
    member_geometry.add_argument("--size-scaling-exponent", type=float, default=0.5)
    member_geometry.add_argument("--minimum-size-scale", type=float, default=0.2)
    member_geometry.add_argument("--source-kind", choices=["gaussian", "beta_model", "softened_nfw"], default="softened_nfw")
    member_geometry.add_argument("--source-amplitude", type=float, default=1.0e-3)
    member_geometry.add_argument("--source-sigma", type=float, default=120.0)
    member_geometry.add_argument("--source-core-radius", type=float, default=60.0)
    member_geometry.add_argument("--source-scale-radius", type=float, default=350.0)
    member_geometry.add_argument("--source-beta", type=float, default=0.8)
    member_geometry.add_argument("--source-outer-slope", type=float, default=3.0)
    member_geometry.add_argument("--screening-mass", type=float, default=0.01)
    member_geometry.add_argument("--density-scale", type=float, default=1.0e-6)
    member_geometry.add_argument("--gravity-scale", type=float, default=1.0)
    member_geometry.add_argument("--exterior-match-tolerance", type=float, default=1.0e-6)
    member_geometry.add_argument("--profile-r-max", type=float, default=8.0e3)
    member_geometry.add_argument("--profile-samples", type=int, default=4096)
    member_geometry.add_argument("--impact-min", type=float, default=30.0)
    member_geometry.add_argument("--impact-max", type=float, default=4.0e3)
    member_geometry.add_argument("--radial-samples", type=int, default=64)
    member_geometry.add_argument("--r-max", type=float, default=8.0e4)
    member_geometry.add_argument("--map-extent", type=float, default=4.0e3)
    member_geometry.add_argument("--grid-size", type=int, default=129)
    member_geometry.add_argument("--radial-bin-width", type=float, default=120.0)
    member_geometry.add_argument("--smooth-sigma-px", type=float, default=3.0)
    member_geometry.add_argument("--axis-samples", type=int, default=257)
    member_geometry.add_argument("--smooth-background-amplitude-fraction", type=float, default=0.0)
    member_geometry.add_argument("--smooth-background-size-multiplier", type=float, default=4.0)
    member_geometry.add_argument("--smooth-background-minimum-size-scale", type=float, default=1.0)
    member_geometry.add_argument("--smooth-background-position-shrink", type=float, default=0.5)
    member_geometry.add_argument("--include-map-arrays", action="store_true")
    member_geometry.set_defaults(func=command_member_geometry_map)

    compare_hff = subparsers.add_parser(
        "compare-hff-member-geometry",
        help="Compare archival HFF residual maps against Einstein member-geometry residual maps.",
    )
    compare_hff.add_argument("--cluster-label", default="")
    compare_hff.add_argument("--member-table", required=True)
    compare_hff.add_argument("--kappa-map", required=True)
    compare_hff.add_argument("--top-members", type=int, default=25)
    compare_hff.add_argument("--mass-key", default="mass_neb")
    compare_hff.add_argument("--member-flag-key", default="is_cluster_member")
    compare_hff.add_argument("--ra-key", default="ra")
    compare_hff.add_argument("--dec-key", default="dec")
    compare_hff.add_argument("--amplitude-scaling-exponent", type=float, default=1.0)
    compare_hff.add_argument("--size-scaling-exponent", type=float, default=0.5)
    compare_hff.add_argument("--minimum-size-scale", type=float, default=0.2)
    compare_hff.add_argument("--source-kind", choices=["gaussian", "beta_model", "softened_nfw"], default="softened_nfw")
    compare_hff.add_argument("--source-amplitude", type=float, default=1.0e-3)
    compare_hff.add_argument("--source-sigma", type=float, default=120.0)
    compare_hff.add_argument("--source-core-radius", type=float, default=60.0)
    compare_hff.add_argument("--source-scale-radius", type=float, default=350.0)
    compare_hff.add_argument("--source-beta", type=float, default=0.8)
    compare_hff.add_argument("--source-outer-slope", type=float, default=3.0)
    compare_hff.add_argument("--screening-mass", type=float, default=0.01)
    compare_hff.add_argument("--density-scale", type=float, default=1.0e-6)
    compare_hff.add_argument("--gravity-scale", type=float, default=1.0)
    compare_hff.add_argument("--exterior-match-tolerance", type=float, default=1.0e-6)
    compare_hff.add_argument("--profile-r-max", type=float, default=8.0e3)
    compare_hff.add_argument("--profile-samples", type=int, default=4096)
    compare_hff.add_argument("--impact-min", type=float, default=30.0)
    compare_hff.add_argument("--impact-max", type=float, default=4.0e3)
    compare_hff.add_argument("--radial-samples", type=int, default=64)
    compare_hff.add_argument("--r-max", type=float, default=8.0e4)
    compare_hff.add_argument("--map-extent", type=float, default=0.0)
    compare_hff.add_argument("--auto-extent-padding-arcsec", type=float, default=80.0)
    compare_hff.add_argument("--minimum-extent-arcsec", type=float, default=120.0)
    compare_hff.add_argument("--grid-size", type=int, default=513)
    compare_hff.add_argument("--radial-bin-width", type=float, default=120.0)
    compare_hff.add_argument("--smooth-sigma-px", type=float, default=3.0)
    compare_hff.add_argument("--axis-samples", type=int, default=257)
    compare_hff.add_argument("--smooth-background-amplitude-fraction", type=float, default=0.0)
    compare_hff.add_argument("--smooth-background-size-multiplier", type=float, default=4.0)
    compare_hff.add_argument("--smooth-background-minimum-size-scale", type=float, default=1.0)
    compare_hff.add_argument("--smooth-background-position-shrink", type=float, default=0.5)
    compare_hff.add_argument("--residual-mode", choices=["gaussian_highpass", "radial_median", "radial_median_bandpass"], default="radial_median_bandpass")
    compare_hff.add_argument("--archive-smooth-sigma-px", type=float, default=18.0)
    compare_hff.add_argument("--archive-radial-bin-arcsec", type=float, default=8.0)
    compare_hff.add_argument("--inner-radius-arcsec", type=float, default=2.5)
    compare_hff.add_argument("--ring-inner-arcsec", type=float, default=5.0)
    compare_hff.add_argument("--ring-outer-arcsec", type=float, default=9.0)
    compare_hff.add_argument("--alignment-max-shift-arcsec", type=float, default=0.0)
    compare_hff.add_argument("--alignment-shift-steps", type=int, default=5)
    compare_hff.add_argument("--alignment-max-rotation-deg", type=float, default=0.0)
    compare_hff.add_argument("--alignment-rotation-steps", type=int, default=5)
    compare_hff.add_argument("--include-map-arrays", action="store_true")
    compare_hff.set_defaults(func=command_compare_hff_member_geometry)

    return parser


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
