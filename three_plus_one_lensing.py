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
import json
import math
import warnings
from dataclasses import dataclass
from typing import Callable

import numpy as np
from scipy.integrate import IntegrationWarning, cumulative_trapezoid, quad


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


def reverse_cumulative_trapezoid(values: np.ndarray, grid: np.ndarray) -> np.ndarray:
    result = np.zeros_like(values)
    if values.size > 1:
        shell_steps = 0.5 * (values[1:] + values[:-1]) * np.diff(grid)
        result[:-1] = np.cumsum(shell_steps[::-1])[::-1]
    return result


def solve_screened_poisson_profile(
    source_amplitude: float,
    source_sigma: float,
    screening_mass: float,
    profile_r_max: float,
    profile_samples: int,
) -> SphericalMeteringProfile:
    if profile_samples < 5:
        raise ValueError("profile_samples must be at least 5.")

    radii = np.linspace(0.0, profile_r_max, profile_samples, dtype=float)
    step = float(radii[1] - radii[0])
    source = gaussian_source(radii, amplitude=source_amplitude, sigma=source_sigma)

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
            "source_amplitude": float(source_amplitude),
            "source_sigma": float(source_sigma),
            "profile_r_max": float(profile_r_max),
            "profile_samples": float(profile_samples),
            "max_mu": float(np.max(mu)),
            "max_source": float(np.max(source)),
            "max_field_equation_residual": float(np.max(np.abs(residual))),
        },
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
    exterior_boundary_phi = 0.5 * math.log(max(1.0 - (2.0 * mass_function[-1] / max(radii[-1], 1.0e-10)), 1.0e-12))
    phi = phi - phi[-1] + exterior_boundary_phi

    lapse_values = np.maximum(np.exp(2.0 * phi), 1.0e-8)
    radial_values = 1.0 / np.maximum(1.0 - compactness, 1.0e-8)
    radial_values[0] = radial_values[1]

    r_min = max(float(radii[1]), 1.0e-6)

    def lapse_A(r: float) -> float:
        return float(np.interp(r, radii, lapse_values, left=lapse_values[1], right=1.0))

    def radial_B(r: float) -> float:
        return float(np.interp(r, radii, radial_values, left=radial_values[1], right=1.0))

    metadata = dict(profile.metadata)
    metadata.update(
        {
            "density_scale": float(density_scale),
            "gravity_scale": float(gravity_scale),
            "max_energy_density": float(np.max(rho)),
            "max_radial_pressure": float(np.max(radial_pressure)),
            "max_compactness_2GM_over_r": float(np.max(compactness)),
            "minimum_photon_lapse": float(np.min(np.sqrt(lapse_values))),
            "gravity_scaled_mass": float(mass_function[-1]),
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

    return parser


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
