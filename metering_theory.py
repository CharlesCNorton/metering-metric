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


@dataclass(frozen=True)
class SourceDomainState:
    activity_density: float
    bridge_conversion: float
    static_source: float


@dataclass(frozen=True)
class BackreactionDomainState:
    activity_density: float
    bridge_conversion: float
    static_source: float
    density_scale: float
    gravity_scale: float
    effective_scale: float


@dataclass(frozen=True)
class MassiveModeObservables:
    rest_mass: float
    conserved_energy: float
    turning_lapse: float
    minimum_lapse: float
    turning_point: float | None
    minimum_group_velocity: float
    mean_group_velocity: float
    crossing_time_s: float
    delay_s: float
    wkb_action: float


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


def activity_density_scalar(number_density: float | np.ndarray, decoherence_rate: float | np.ndarray) -> float | np.ndarray:
    return np.asarray(number_density) * np.asarray(decoherence_rate)


def occupancy_bridge_conversion(eta_j: float, tau_p_s: float, l_perp_m: float) -> float:
    return OccupancyBridge(eta_j=eta_j, tau_p_s=tau_p_s, l_perp_m=l_perp_m).conversion_factor()


def bridge_efficiency_for_conversion(required_conversion: float, tau_p_s: float, l_perp_m: float) -> float:
    if required_conversion <= 0.0:
        raise ValueError("required_conversion must be positive.")
    if tau_p_s <= 0.0:
        raise ValueError("tau_p_s must be positive.")
    if l_perp_m <= 0.0:
        raise ValueError("l_perp_m must be positive.")
    return float(required_conversion * (l_perp_m * l_perp_m) / tau_p_s)


def bridge_transverse_scale_for_conversion(required_conversion: float, eta_j: float, tau_p_s: float) -> float:
    if required_conversion <= 0.0:
        raise ValueError("required_conversion must be positive.")
    if eta_j <= 0.0:
        raise ValueError("eta_j must be positive.")
    if tau_p_s <= 0.0:
        raise ValueError("tau_p_s must be positive.")
    return float(math.sqrt(eta_j * tau_p_s / required_conversion))


def bridge_persistence_for_conversion(required_conversion: float, eta_j: float, l_perp_m: float) -> float:
    if required_conversion <= 0.0:
        raise ValueError("required_conversion must be positive.")
    if eta_j <= 0.0:
        raise ValueError("eta_j must be positive.")
    if l_perp_m <= 0.0:
        raise ValueError("l_perp_m must be positive.")
    return float(required_conversion * (l_perp_m * l_perp_m) / eta_j)


def effective_backreaction_scale(source_amplitude: float, density_scale: float, gravity_scale: float) -> float:
    if source_amplitude <= 0.0:
        raise ValueError("source_amplitude must be positive.")
    if density_scale <= 0.0:
        raise ValueError("density_scale must be positive.")
    if gravity_scale <= 0.0:
        raise ValueError("gravity_scale must be positive.")
    return float(gravity_scale * density_scale * source_amplitude * source_amplitude)


def source_domain_state(activity_density: float, bridge_conversion: float) -> SourceDomainState:
    static_source = static_source_amplitude_from_activity(
        activity_density=activity_density,
        bridge_conversion=bridge_conversion,
    )
    return SourceDomainState(
        activity_density=float(activity_density),
        bridge_conversion=float(bridge_conversion),
        static_source=float(static_source),
    )


def backreaction_domain_state(
    activity_density: float,
    bridge_conversion: float,
    density_scale: float,
    gravity_scale: float,
) -> BackreactionDomainState:
    source_state = source_domain_state(
        activity_density=activity_density,
        bridge_conversion=bridge_conversion,
    )
    effective_scale = effective_backreaction_scale(
        source_amplitude=source_state.static_source,
        density_scale=density_scale,
        gravity_scale=gravity_scale,
    )
    return BackreactionDomainState(
        activity_density=source_state.activity_density,
        bridge_conversion=source_state.bridge_conversion,
        static_source=source_state.static_source,
        density_scale=float(density_scale),
        gravity_scale=float(gravity_scale),
        effective_scale=float(effective_scale),
    )


def static_source_amplitude_from_activity(activity_density: float, bridge_conversion: float) -> float:
    if activity_density <= 0.0:
        raise ValueError("activity_density must be positive.")
    if bridge_conversion <= 0.0:
        raise ValueError("bridge_conversion must be positive.")
    return float(bridge_conversion * activity_density)


def backreaction_scale_from_activity(
    activity_density: float,
    bridge_conversion: float,
    density_scale: float,
    gravity_scale: float,
) -> float:
    source_amplitude = static_source_amplitude_from_activity(
        activity_density=activity_density,
        bridge_conversion=bridge_conversion,
    )
    return effective_backreaction_scale(
        source_amplitude=source_amplitude,
        density_scale=density_scale,
        gravity_scale=gravity_scale,
    )


def required_bridge_conversion_from_backreaction_scale(
    effective_scale: float,
    activity_density: float,
    density_scale: float,
    gravity_scale: float,
) -> float:
    if effective_scale <= 0.0:
        raise ValueError("effective_scale must be positive.")
    if activity_density <= 0.0:
        raise ValueError("activity_density must be positive.")
    if density_scale <= 0.0:
        raise ValueError("density_scale must be positive.")
    if gravity_scale <= 0.0:
        raise ValueError("gravity_scale must be positive.")
    return float(math.sqrt(effective_scale / (gravity_scale * density_scale)) / activity_density)


def required_activity_density_from_backreaction_scale(
    effective_scale: float,
    bridge_conversion: float,
    density_scale: float,
    gravity_scale: float,
) -> float:
    if effective_scale <= 0.0:
        raise ValueError("effective_scale must be positive.")
    if bridge_conversion <= 0.0:
        raise ValueError("bridge_conversion must be positive.")
    if density_scale <= 0.0:
        raise ValueError("density_scale must be positive.")
    if gravity_scale <= 0.0:
        raise ValueError("gravity_scale must be positive.")
    return float(math.sqrt(effective_scale / (gravity_scale * density_scale)) / bridge_conversion)


def static_source_from_activity(activity_density: float | np.ndarray, bridge: OccupancyBridge) -> float | np.ndarray:
    return bridge.conversion_factor() * np.asarray(activity_density)


def source_ratio_from_activity(
    activity_density_numerator: float,
    activity_density_denominator: float,
) -> float:
    if activity_density_denominator == 0.0:
        raise ValueError("activity_density_denominator must be nonzero.")
    return float(activity_density_numerator / activity_density_denominator)


def bridge_matched_source_amplitude(
    reference_source_amplitude: float,
    activity_density: float,
    reference_activity_density: float,
) -> float:
    if reference_source_amplitude <= 0.0:
        raise ValueError("reference_source_amplitude must be positive.")
    ratio = source_ratio_from_activity(
        activity_density_numerator=activity_density,
        activity_density_denominator=reference_activity_density,
    )
    return float(reference_source_amplitude * ratio)


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


def same_lapse_profile_crossing_time_seconds(
    grid_m: np.ndarray,
    lapse: np.ndarray,
    conserved_energy: float,
    rest_mass: float = 1.0,
) -> float:
    grid_m = np.asarray(grid_m, dtype=float)
    lapse = np.asarray(lapse, dtype=float)
    if grid_m.ndim != 1 or lapse.ndim != 1 or grid_m.size != lapse.size:
        raise ValueError("grid_m and lapse must be 1D arrays of equal length.")
    velocity = same_lapse_group_velocity(
        lapse=lapse,
        conserved_energy=conserved_energy,
        rest_mass=rest_mass,
    )
    if np.any(velocity <= 1.0e-15):
        return math.inf
    return float(np.trapezoid(1.0 / (C * velocity), grid_m))


def same_lapse_profile_delay_seconds(
    grid_m: np.ndarray,
    lapse: np.ndarray,
    conserved_energy: float,
    rest_mass: float = 1.0,
) -> float:
    crossing_time = same_lapse_profile_crossing_time_seconds(
        grid_m=grid_m,
        lapse=lapse,
        conserved_energy=conserved_energy,
        rest_mass=rest_mass,
    )
    if not math.isfinite(crossing_time):
        return math.inf
    grid_m = np.asarray(grid_m, dtype=float)
    return float(crossing_time - (grid_m[-1] - grid_m[0]) / C)


def same_lapse_profile_wkb_action(
    grid: np.ndarray,
    lapse: np.ndarray,
    conserved_energy: float,
    rest_mass: float = 1.0,
) -> float:
    density = same_lapse_wkb_action_density(
        lapse=lapse,
        conserved_energy=conserved_energy,
        rest_mass=rest_mass,
    )
    return integrate_positive_density(grid=np.asarray(grid, dtype=float), density=density)


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


def same_lapse_profile_observables(
    grid_m: np.ndarray,
    lapse: np.ndarray,
    conserved_energy: float,
    rest_mass: float = 1.0,
) -> MassiveModeObservables:
    grid_m = np.asarray(grid_m, dtype=float)
    lapse = np.asarray(lapse, dtype=float)
    if grid_m.ndim != 1 or lapse.ndim != 1 or grid_m.size != lapse.size:
        raise ValueError("grid_m and lapse must be 1D arrays of equal length.")
    velocity = same_lapse_group_velocity(
        lapse=lapse,
        conserved_energy=conserved_energy,
        rest_mass=rest_mass,
    )
    finite_velocity = velocity[np.isfinite(velocity)]
    minimum_velocity = float(np.min(finite_velocity)) if finite_velocity.size else float("nan")
    mean_velocity = float(np.mean(finite_velocity)) if finite_velocity.size else float("nan")
    return MassiveModeObservables(
        rest_mass=float(rest_mass),
        conserved_energy=float(conserved_energy),
        turning_lapse=float(same_lapse_turning_lapse(conserved_energy=conserved_energy, rest_mass=rest_mass)),
        minimum_lapse=float(np.min(lapse)),
        turning_point=first_turning_point(
            grid=grid_m,
            lapse=lapse,
            conserved_energy=conserved_energy,
            rest_mass=rest_mass,
        ),
        minimum_group_velocity=minimum_velocity,
        mean_group_velocity=mean_velocity,
        crossing_time_s=same_lapse_profile_crossing_time_seconds(
            grid_m=grid_m,
            lapse=lapse,
            conserved_energy=conserved_energy,
            rest_mass=rest_mass,
        ),
        delay_s=same_lapse_profile_delay_seconds(
            grid_m=grid_m,
            lapse=lapse,
            conserved_energy=conserved_energy,
            rest_mass=rest_mass,
        ),
        wkb_action=same_lapse_profile_wkb_action(
            grid=grid_m,
            lapse=lapse,
            conserved_energy=conserved_energy,
            rest_mass=rest_mass,
        ),
    )


def dirichlet_schrodinger_matrix(grid: np.ndarray, potential: np.ndarray) -> np.ndarray:
    grid = np.asarray(grid, dtype=float)
    potential = np.asarray(potential, dtype=float)
    if grid.ndim != 1 or potential.ndim != 1 or grid.size != potential.size:
        raise ValueError("grid and potential must be 1D arrays of equal length.")
    if grid.size < 3:
        raise ValueError("grid must have at least three points.")
    spacing = np.diff(grid)
    if not np.allclose(spacing, spacing[0], rtol=1.0e-9, atol=1.0e-12):
        raise ValueError("grid must be uniformly spaced.")
    dx = float(spacing[0])
    interior_potential = potential[1:-1]
    diagonal = (2.0 / (dx * dx)) + interior_potential
    off_diagonal = np.full(interior_potential.size - 1, -1.0 / (dx * dx), dtype=float)
    matrix = np.diag(diagonal)
    if off_diagonal.size > 0:
        matrix += np.diag(off_diagonal, 1) + np.diag(off_diagonal, -1)
    return matrix


def lowest_dirichlet_eigenvalues(
    grid: np.ndarray,
    potential: np.ndarray,
    count: int,
) -> np.ndarray:
    if count <= 0:
        raise ValueError("count must be positive.")
    matrix = dirichlet_schrodinger_matrix(grid=grid, potential=potential)
    eigenvalues = np.linalg.eigvalsh(matrix)
    return np.asarray(eigenvalues[:count], dtype=float)
