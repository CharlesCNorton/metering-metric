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
from dataclasses import dataclass
from typing import Callable

import numpy as np
from scipy.integrate import quad


MetricFn = Callable[[float], float]


@dataclass(frozen=True)
class StaticSphericalPhotonMetric:
    name: str
    lapse_A: MetricFn
    radial_B: MetricFn
    r_min: float
    asymptotic_lapse: float = 1.0


def tanh_lapse(mu: float, alpha: float, epsilon: float = 0.0) -> float:
    return float(epsilon + (1.0 - epsilon) * math.tanh(alpha * mu))


def gaussian_mu(r: float, mu0: float, sigma: float) -> float:
    return float(mu0 * math.exp(-(r * r) / (2.0 * sigma * sigma)))


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


def photon_deflection_angle(
    metric: StaticSphericalPhotonMetric,
    impact_parameter: float,
    r_max: float,
) -> dict:
    turning_radius = find_turning_radius(metric, impact_parameter=impact_parameter, r_max=r_max)
    metric_integral, metric_error = quad(
        lambda r: deflection_integrand(metric, r, impact_parameter),
        turning_radius,
        r_max,
        points=[turning_radius],
        limit=400,
        epsabs=1.0e-10,
        epsrel=1.0e-10,
    )
    flat_integral, flat_error = quad(
        lambda r: deflection_integrand(flat_decoupled_metric(), r, impact_parameter),
        impact_parameter,
        r_max,
        points=[impact_parameter],
        limit=400,
        epsabs=1.0e-10,
        epsrel=1.0e-10,
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
    metric_time, metric_error = quad(
        lambda r: coordinate_time_integrand(metric, r, impact_parameter),
        turning_radius,
        r_observer,
        points=[turning_radius],
        limit=400,
        epsabs=1.0e-10,
        epsrel=1.0e-10,
    )
    flat_time, flat_error = quad(
        lambda r: coordinate_time_integrand(flat_decoupled_metric(), r, impact_parameter),
        impact_parameter,
        r_observer,
        points=[impact_parameter],
        limit=400,
        epsabs=1.0e-10,
        epsrel=1.0e-10,
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


def scenario_metric(args: argparse.Namespace) -> StaticSphericalPhotonMetric:
    if args.scenario == "flat_decoupled":
        return flat_decoupled_metric()
    if args.scenario == "schwarzschild":
        return schwarzschild_metric(mass=args.mass)
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
    print(json.dumps(result, indent=2))


def command_benchmark(args: argparse.Namespace) -> None:
    result = weak_field_schwarzschild_benchmark(
        mass=args.mass,
        impact_parameter=args.impact_parameter,
        r_max=args.r_max,
    )
    print(json.dumps(result, indent=2))


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Evaluate 3+1D photon observables for static spherical metrics.",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    evaluate = subparsers.add_parser("evaluate", help="Evaluate deflection, time delay, and redshift for a metric.")
    evaluate.add_argument("--scenario", choices=["flat_decoupled", "schwarzschild", "toy_photon_coupled_metering"], required=True)
    evaluate.add_argument("--impact-parameter", type=float, default=50.0)
    evaluate.add_argument("--r-max", type=float, default=5.0e4)
    evaluate.add_argument("--r-emit", type=float, default=10.0)
    evaluate.add_argument("--r-obs", type=float, default=1.0e6)
    evaluate.add_argument("--mass", type=float, default=1.0)
    evaluate.add_argument("--mu0", type=float, default=1.0)
    evaluate.add_argument("--sigma", type=float, default=100.0)
    evaluate.add_argument("--alpha", type=float, default=5.0)
    evaluate.add_argument("--epsilon", type=float, default=0.05)
    evaluate.set_defaults(func=command_evaluate)

    benchmark = subparsers.add_parser("benchmark-schwarzschild", help="Check the null-geodesic integral against 4M/b.")
    benchmark.add_argument("--mass", type=float, default=1.0)
    benchmark.add_argument("--impact-parameter", type=float, default=1000.0)
    benchmark.add_argument("--r-max", type=float, default=1.0e6)
    benchmark.set_defaults(func=command_benchmark)

    return parser


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
