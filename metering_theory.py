from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np


HBAR = 1.054_571_817e-34
C = 299_792_458.0
G_NEWTON = 6.674_30e-11
PLANCK_TIME = math.sqrt(HBAR * G_NEWTON / (C**5))


@dataclass(frozen=True)
class AlphaRegimes:
    lattice_fit: float
    continuum_upper: float
    laboratory_benchmark_upper: float
    physical_si_estimate: float


@dataclass(frozen=True)
class OccupancyBridge:
    eta_j: float
    tau_p_s: float
    l_perp_m: float

    def conversion_factor(self) -> float:
        if self.tau_p_s <= 0.0:
            raise ValueError("tau_p_s must be positive.")
        if self.l_perp_m <= 0.0:
            raise ValueError("l_perp_m must be positive.")
        return float(self.eta_j * self.tau_p_s / (self.l_perp_m * self.l_perp_m))


def tanh_lapse(mu: float | np.ndarray, alpha: float, epsilon: float = 0.0) -> float | np.ndarray:
    if epsilon < 0.0 or epsilon >= 1.0:
        raise ValueError("epsilon must lie in [0, 1).")
    result = epsilon + (1.0 - epsilon) * np.tanh(alpha * np.asarray(mu))
    if np.ndim(result) == 0:
        return float(result)
    return result


def gaussian_metering_density(grid: np.ndarray, mu0: float, sigma: float) -> np.ndarray:
    if sigma <= 0.0:
        raise ValueError("sigma must be positive.")
    return mu0 * np.exp(-(grid * grid) / (2.0 * sigma * sigma))


def occupancy_bridge_conversion(eta_j: float, tau_p_s: float, l_perp_m: float) -> float:
    return OccupancyBridge(eta_j=eta_j, tau_p_s=tau_p_s, l_perp_m=l_perp_m).conversion_factor()


def static_source_from_activity(activity_density: float | np.ndarray, bridge: OccupancyBridge) -> float | np.ndarray:
    return bridge.conversion_factor() * np.asarray(activity_density)


def uniform_screened_mu_from_source(static_source: float | np.ndarray, screening_mass: float) -> float | np.ndarray:
    if screening_mass <= 0.0:
        raise ValueError("screening_mass must be positive.")
    return np.asarray(static_source) / (screening_mass * screening_mass)


def same_lapse_effective_potential(lapse: np.ndarray, rest_mass: float = 1.0) -> np.ndarray:
    lapse_safe = np.maximum(np.asarray(lapse, dtype=float), 1.0e-18)
    return (rest_mass / lapse_safe) ** 2


def same_lapse_wave_number(lapse: np.ndarray, conserved_energy: float, rest_mass: float = 1.0) -> np.ndarray:
    potential = same_lapse_effective_potential(lapse=lapse, rest_mass=rest_mass)
    return np.sqrt(np.maximum(conserved_energy * conserved_energy - potential, 0.0))


def same_lapse_group_velocity(lapse: np.ndarray, conserved_energy: float, rest_mass: float = 1.0) -> np.ndarray:
    if conserved_energy <= 0.0:
        raise ValueError("conserved_energy must be positive.")
    wave_number = same_lapse_wave_number(lapse=lapse, conserved_energy=conserved_energy, rest_mass=rest_mass)
    return wave_number / conserved_energy


def same_lapse_delay_density(lapse: np.ndarray, conserved_energy: float, rest_mass: float = 1.0) -> np.ndarray:
    velocity = same_lapse_group_velocity(lapse=lapse, conserved_energy=conserved_energy, rest_mass=rest_mass)
    result = np.full_like(velocity, np.inf, dtype=float)
    mask = velocity > 1.0e-15
    result[mask] = (1.0 / velocity[mask]) - 1.0
    return result


def same_lapse_tortoise_density(lapse: np.ndarray) -> np.ndarray:
    lapse_safe = np.maximum(np.asarray(lapse, dtype=float), 1.0e-18)
    return 1.0 / lapse_safe


def same_lapse_turning_lapse(conserved_energy: float, rest_mass: float = 1.0) -> float:
    if conserved_energy <= 0.0:
        raise ValueError("conserved_energy must be positive.")
    return rest_mass / conserved_energy


def same_lapse_wkb_action_density(lapse: np.ndarray, conserved_energy: float, rest_mass: float = 1.0) -> np.ndarray:
    potential = same_lapse_effective_potential(lapse=lapse, rest_mass=rest_mass)
    return np.sqrt(np.maximum(potential - conserved_energy * conserved_energy, 0.0))


def same_lapse_uniform_crossing_time_seconds(
    path_length_m: float,
    lapse: float,
    conserved_energy: float,
    rest_mass: float = 1.0,
) -> float:
    if path_length_m < 0.0:
        raise ValueError("path_length_m must be nonnegative.")
    velocity = same_lapse_group_velocity(
        lapse=np.asarray([lapse], dtype=float),
        conserved_energy=conserved_energy,
        rest_mass=rest_mass,
    )[0]
    if velocity <= 1.0e-15:
        return math.inf
    return float(path_length_m / (C * velocity))


def same_lapse_uniform_delay_seconds(
    path_length_m: float,
    lapse: float,
    conserved_energy: float,
    rest_mass: float = 1.0,
) -> float:
    crossing_time = same_lapse_uniform_crossing_time_seconds(
        path_length_m=path_length_m,
        lapse=lapse,
        conserved_energy=conserved_energy,
        rest_mass=rest_mass,
    )
    if not math.isfinite(crossing_time):
        return math.inf
    return float(crossing_time - path_length_m / C)


def integrate_positive_density(grid: np.ndarray, density: np.ndarray) -> float:
    grid = np.asarray(grid, dtype=float)
    density = np.asarray(density, dtype=float)
    if grid.ndim != 1 or density.ndim != 1 or grid.size != density.size:
        raise ValueError("grid and density must be 1D arrays of equal length.")
    return float(np.trapezoid(np.maximum(density, 0.0), grid))


def first_turning_point(grid: np.ndarray, lapse: np.ndarray, conserved_energy: float, rest_mass: float = 1.0) -> float | None:
    grid = np.asarray(grid, dtype=float)
    lapse = np.asarray(lapse, dtype=float)
    threshold = same_lapse_turning_lapse(conserved_energy=conserved_energy, rest_mass=rest_mass)
    allowed = lapse >= threshold
    if np.all(allowed):
        return None
    first_forbidden = int(np.argmax(~allowed))
    if first_forbidden == 0:
        return float(grid[0])
    x0 = grid[first_forbidden - 1]
    x1 = grid[first_forbidden]
    n0 = lapse[first_forbidden - 1]
    n1 = lapse[first_forbidden]
    if abs(n1 - n0) < 1.0e-18:
        return float(x1)
    weight = (threshold - n0) / (n1 - n0)
    return float(x0 + weight * (x1 - x0))
