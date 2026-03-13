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
class TheoryNativeBridgeClosure:
    screening_length_m: float
    eta_j: float
    tau_planck_multiples: float = 1.0
    transverse_screening_multiples: float = 1.0
    expansion_rate_s_inv: float = 0.0

    def tau_p_s(self) -> float:
        if self.tau_planck_multiples <= 0.0:
            raise ValueError("tau_planck_multiples must be positive.")
        return float(self.tau_planck_multiples * PLANCK_TIME)

    def l_perp_m(self) -> float:
        if self.screening_length_m <= 0.0:
            raise ValueError("screening_length_m must be positive.")
        if self.transverse_screening_multiples <= 0.0:
            raise ValueError("transverse_screening_multiples must be positive.")
        return float(self.transverse_screening_multiples * self.screening_length_m)

    def dimensionless_closure_ratio(self) -> float:
        if self.eta_j <= 0.0:
            raise ValueError("eta_j must be positive.")
        if self.tau_planck_multiples <= 0.0:
            raise ValueError("tau_planck_multiples must be positive.")
        if self.transverse_screening_multiples <= 0.0:
            raise ValueError("transverse_screening_multiples must be positive.")
        return float(self.eta_j * self.tau_planck_multiples / (self.transverse_screening_multiples ** 2))

    def bridge_conversion(self) -> float:
        return relaxed_bridge_conversion(
            eta_j=self.eta_j,
            tau_p_s=self.tau_p_s(),
            l_perp_m=self.l_perp_m(),
            expansion_rate_s_inv=self.expansion_rate_s_inv,
        )


@dataclass(frozen=True)
class OccupancyRelaxationState:
    activity_density: float
    expansion_rate_s_inv: float
    tau_p_s: float
    occupancy_density: float
    eta_j: float
    l_perp_m: float
    bridge_conversion: float
    static_source: float


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


@dataclass(frozen=True)
class SpeciesProfileObservables:
    species_name: str
    rest_mass: float
    conserved_energy: float
    turning_lapse: float
    required_activity_density: float
    observables: MassiveModeObservables


@dataclass(frozen=True)
class PhysicalSpeciesKinematics:
    species_name: str
    rest_mass_eV: float
    conserved_energy_eV: float

    def mass_energy_ratio(self) -> float:
        if self.rest_mass_eV < 0.0:
            raise ValueError("rest_mass_eV must be nonnegative.")
        if self.conserved_energy_eV <= 0.0:
            raise ValueError("conserved_energy_eV must be positive.")
        return float(self.rest_mass_eV / self.conserved_energy_eV)


@dataclass(frozen=True)
class TheoryNativeActivityDomain:
    domain_name: str
    screening_length_m: float
    activity_density: float
    path_length_m: float
    alpha: float
    closure_ratio: float
    epsilon: float = 0.0
    expansion_rate_s_inv: float = 0.0
    tau_planck_multiples: float = 1.0


@dataclass(frozen=True)
class DomainSpeciesPrediction:
    domain_name: str
    species_name: str
    bridge_conversion: float
    closure_ratio: float
    uniform_lapse: float
    mass_energy_ratio: float
    required_activity_density: float
    excluded: bool
    observables: MassiveModeObservables


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


def occupancy_relaxation_factor(tau_p_s: float, expansion_rate_s_inv: float = 0.0) -> float:
    if tau_p_s <= 0.0:
        raise ValueError("tau_p_s must be positive.")
    denominator = 1.0 + (expansion_rate_s_inv * tau_p_s)
    if denominator <= 0.0:
        raise ValueError("1 + expansion_rate_s_inv * tau_p_s must be positive.")
    return float(1.0 / denominator)


def steady_state_occupancy_density(
    activity_density: float,
    tau_p_s: float,
    expansion_rate_s_inv: float = 0.0,
) -> float:
    if activity_density <= 0.0:
        raise ValueError("activity_density must be positive.")
    relaxation_factor = occupancy_relaxation_factor(
        tau_p_s=tau_p_s,
        expansion_rate_s_inv=expansion_rate_s_inv,
    )
    return float(activity_density * tau_p_s * relaxation_factor)


def relaxed_bridge_conversion(
    eta_j: float,
    tau_p_s: float,
    l_perp_m: float,
    expansion_rate_s_inv: float = 0.0,
) -> float:
    if eta_j <= 0.0:
        raise ValueError("eta_j must be positive.")
    if l_perp_m <= 0.0:
        raise ValueError("l_perp_m must be positive.")
    relaxation_factor = occupancy_relaxation_factor(
        tau_p_s=tau_p_s,
        expansion_rate_s_inv=expansion_rate_s_inv,
    )
    return float(eta_j * tau_p_s * relaxation_factor / (l_perp_m * l_perp_m))


def theory_native_bridge_conversion(
    screening_length_m: float,
    eta_j: float,
    tau_planck_multiples: float = 1.0,
    transverse_screening_multiples: float = 1.0,
    expansion_rate_s_inv: float = 0.0,
) -> float:
    closure = TheoryNativeBridgeClosure(
        screening_length_m=screening_length_m,
        eta_j=eta_j,
        tau_planck_multiples=tau_planck_multiples,
        transverse_screening_multiples=transverse_screening_multiples,
        expansion_rate_s_inv=expansion_rate_s_inv,
    )
    return closure.bridge_conversion()


def required_theory_native_closure_ratio(
    required_conversion: float,
    screening_length_m: float,
    expansion_rate_s_inv: float = 0.0,
    tau_planck_multiples: float = 1.0,
) -> float:
    if required_conversion <= 0.0:
        raise ValueError("required_conversion must be positive.")
    if screening_length_m <= 0.0:
        raise ValueError("screening_length_m must be positive.")
    if tau_planck_multiples <= 0.0:
        raise ValueError("tau_planck_multiples must be positive.")
    relaxation_factor = occupancy_relaxation_factor(
        tau_p_s=tau_planck_multiples * PLANCK_TIME,
        expansion_rate_s_inv=expansion_rate_s_inv,
    )
    return float(required_conversion * (screening_length_m * screening_length_m) / (PLANCK_TIME * relaxation_factor))


def bridge_conversion_from_theory_native_closure_ratio(
    closure_ratio: float,
    screening_length_m: float,
    expansion_rate_s_inv: float = 0.0,
    tau_planck_multiples: float = 1.0,
) -> float:
    if closure_ratio <= 0.0:
        raise ValueError("closure_ratio must be positive.")
    if screening_length_m <= 0.0:
        raise ValueError("screening_length_m must be positive.")
    if tau_planck_multiples <= 0.0:
        raise ValueError("tau_planck_multiples must be positive.")
    relaxation_factor = occupancy_relaxation_factor(
        tau_p_s=tau_planck_multiples * PLANCK_TIME,
        expansion_rate_s_inv=expansion_rate_s_inv,
    )
    return float(closure_ratio * PLANCK_TIME * relaxation_factor / (screening_length_m * screening_length_m))


def required_screening_cell_occupancy_fraction(
    required_conversion: float,
    screening_length_m: float,
    expansion_rate_s_inv: float = 0.0,
    tau_planck_multiples: float = 1.0,
) -> float:
    return required_theory_native_closure_ratio(
        required_conversion=required_conversion,
        screening_length_m=screening_length_m,
        expansion_rate_s_inv=expansion_rate_s_inv,
        tau_planck_multiples=tau_planck_multiples,
    )


def theory_native_transverse_multiple_for_unit_ratio(
    required_conversion: float,
    screening_length_m: float,
    expansion_rate_s_inv: float = 0.0,
    tau_planck_multiples: float = 1.0,
) -> float:
    ratio = required_theory_native_closure_ratio(
        required_conversion=required_conversion,
        screening_length_m=screening_length_m,
        expansion_rate_s_inv=expansion_rate_s_inv,
        tau_planck_multiples=tau_planck_multiples,
    )
    return float(math.sqrt(max(1.0 / ratio, 0.0)))


def required_eta_for_theory_native_closure(
    required_conversion: float,
    screening_length_m: float,
    tau_planck_multiples: float = 1.0,
    transverse_screening_multiples: float = 1.0,
    expansion_rate_s_inv: float = 0.0,
) -> float:
    ratio = required_theory_native_closure_ratio(
        required_conversion=required_conversion,
        screening_length_m=screening_length_m,
        expansion_rate_s_inv=expansion_rate_s_inv,
        tau_planck_multiples=tau_planck_multiples,
    )
    if tau_planck_multiples <= 0.0:
        raise ValueError("tau_planck_multiples must be positive.")
    if transverse_screening_multiples <= 0.0:
        raise ValueError("transverse_screening_multiples must be positive.")
    return float(ratio * (transverse_screening_multiples ** 2) / tau_planck_multiples)


def occupancy_relaxation_state(
    activity_density: float,
    eta_j: float,
    tau_p_s: float,
    l_perp_m: float,
    expansion_rate_s_inv: float = 0.0,
) -> OccupancyRelaxationState:
    occupancy_density = steady_state_occupancy_density(
        activity_density=activity_density,
        tau_p_s=tau_p_s,
        expansion_rate_s_inv=expansion_rate_s_inv,
    )
    bridge_conversion = relaxed_bridge_conversion(
        eta_j=eta_j,
        tau_p_s=tau_p_s,
        l_perp_m=l_perp_m,
        expansion_rate_s_inv=expansion_rate_s_inv,
    )
    static_source = static_source_amplitude_from_activity(
        activity_density=activity_density,
        bridge_conversion=bridge_conversion,
    )
    return OccupancyRelaxationState(
        activity_density=float(activity_density),
        expansion_rate_s_inv=float(expansion_rate_s_inv),
        tau_p_s=float(tau_p_s),
        occupancy_density=float(occupancy_density),
        eta_j=float(eta_j),
        l_perp_m=float(l_perp_m),
        bridge_conversion=float(bridge_conversion),
        static_source=float(static_source),
    )


def bridge_efficiency_for_conversion(required_conversion: float, tau_p_s: float, l_perp_m: float) -> float:
    if required_conversion <= 0.0:
        raise ValueError("required_conversion must be positive.")
    if tau_p_s <= 0.0:
        raise ValueError("tau_p_s must be positive.")
    if l_perp_m <= 0.0:
        raise ValueError("l_perp_m must be positive.")
    return float(required_conversion * (l_perp_m * l_perp_m) / tau_p_s)


def relaxed_bridge_efficiency_for_conversion(
    required_conversion: float,
    tau_p_s: float,
    l_perp_m: float,
    expansion_rate_s_inv: float = 0.0,
) -> float:
    if required_conversion <= 0.0:
        raise ValueError("required_conversion must be positive.")
    relaxation_factor = occupancy_relaxation_factor(
        tau_p_s=tau_p_s,
        expansion_rate_s_inv=expansion_rate_s_inv,
    )
    return float(required_conversion * (l_perp_m * l_perp_m) / (tau_p_s * relaxation_factor))


def bridge_transverse_scale_for_conversion(required_conversion: float, eta_j: float, tau_p_s: float) -> float:
    if required_conversion <= 0.0:
        raise ValueError("required_conversion must be positive.")
    if eta_j <= 0.0:
        raise ValueError("eta_j must be positive.")
    if tau_p_s <= 0.0:
        raise ValueError("tau_p_s must be positive.")
    return float(math.sqrt(eta_j * tau_p_s / required_conversion))


def relaxed_bridge_transverse_scale_for_conversion(
    required_conversion: float,
    eta_j: float,
    tau_p_s: float,
    expansion_rate_s_inv: float = 0.0,
) -> float:
    if required_conversion <= 0.0:
        raise ValueError("required_conversion must be positive.")
    if eta_j <= 0.0:
        raise ValueError("eta_j must be positive.")
    relaxation_factor = occupancy_relaxation_factor(
        tau_p_s=tau_p_s,
        expansion_rate_s_inv=expansion_rate_s_inv,
    )
    return float(math.sqrt(eta_j * tau_p_s * relaxation_factor / required_conversion))


def bridge_persistence_for_conversion(required_conversion: float, eta_j: float, l_perp_m: float) -> float:
    if required_conversion <= 0.0:
        raise ValueError("required_conversion must be positive.")
    if eta_j <= 0.0:
        raise ValueError("eta_j must be positive.")
    if l_perp_m <= 0.0:
        raise ValueError("l_perp_m must be positive.")
    return float(required_conversion * (l_perp_m * l_perp_m) / eta_j)


def relaxed_bridge_persistence_for_conversion(
    required_conversion: float,
    eta_j: float,
    l_perp_m: float,
    expansion_rate_s_inv: float = 0.0,
) -> float:
    if required_conversion <= 0.0:
        raise ValueError("required_conversion must be positive.")
    if eta_j <= 0.0:
        raise ValueError("eta_j must be positive.")
    if l_perp_m <= 0.0:
        raise ValueError("l_perp_m must be positive.")
    if expansion_rate_s_inv == 0.0:
        return bridge_persistence_for_conversion(
            required_conversion=required_conversion,
            eta_j=eta_j,
            l_perp_m=l_perp_m,
        )
    a = required_conversion / eta_j
    b = required_conversion * expansion_rate_s_inv * (l_perp_m * l_perp_m) / eta_j
    if abs(b) <= 1.0e-30:
        return float(a * l_perp_m * l_perp_m)
    tau_p_s = a * l_perp_m * l_perp_m / max(1.0 - b, 1.0e-30)
    if tau_p_s <= 0.0:
        raise ValueError("No positive persistence time satisfies the relaxed bridge parameters.")
    return float(tau_p_s)


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


def required_activity_density_for_source_amplitude(source_amplitude: float, bridge_conversion: float) -> float:
    if source_amplitude <= 0.0:
        raise ValueError("source_amplitude must be positive.")
    if bridge_conversion <= 0.0:
        raise ValueError("bridge_conversion must be positive.")
    return float(source_amplitude / bridge_conversion)


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


def required_gravity_density_product(
    effective_scale: float,
    source_amplitude: float,
) -> float:
    if effective_scale <= 0.0:
        raise ValueError("effective_scale must be positive.")
    if source_amplitude <= 0.0:
        raise ValueError("source_amplitude must be positive.")
    return float(effective_scale / (source_amplitude * source_amplitude))


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


def required_uniform_mu_for_lapse(target_lapse: float, alpha: float, epsilon: float = 0.0) -> float:
    if alpha <= 0.0:
        raise ValueError("alpha must be positive.")
    if epsilon < 0.0 or epsilon >= 1.0:
        raise ValueError("epsilon must lie in [0, 1).")
    if target_lapse <= epsilon:
        return math.inf
    if target_lapse >= 1.0:
        return math.inf
    reduced = (target_lapse - epsilon) / (1.0 - epsilon)
    if reduced <= 0.0 or reduced >= 1.0:
        return math.inf
    return float(np.arctanh(reduced) / alpha)


def required_uniform_source_for_lapse(
    target_lapse: float,
    alpha: float,
    screening_mass: float,
    epsilon: float = 0.0,
) -> float:
    mu_value = required_uniform_mu_for_lapse(
        target_lapse=target_lapse,
        alpha=alpha,
        epsilon=epsilon,
    )
    return float((screening_mass * screening_mass) * mu_value)


def required_activity_density_for_lapse(
    target_lapse: float,
    alpha: float,
    screening_mass: float,
    bridge_conversion: float,
    epsilon: float = 0.0,
) -> float:
    if bridge_conversion <= 0.0:
        raise ValueError("bridge_conversion must be positive.")
    source_value = required_uniform_source_for_lapse(
        target_lapse=target_lapse,
        alpha=alpha,
        screening_mass=screening_mass,
        epsilon=epsilon,
    )
    if not math.isfinite(source_value):
        return math.inf
    return float(source_value / bridge_conversion)


def required_activity_density_for_turning(
    conserved_energy: float,
    rest_mass: float,
    alpha: float,
    screening_mass: float,
    bridge_conversion: float,
    epsilon: float = 0.0,
) -> float:
    target_lapse = same_lapse_turning_lapse(
        conserved_energy=conserved_energy,
        rest_mass=rest_mass,
    )
    return required_activity_density_for_lapse(
        target_lapse=target_lapse,
        alpha=alpha,
        screening_mass=screening_mass,
        bridge_conversion=bridge_conversion,
        epsilon=epsilon,
    )


def uniform_lapse_from_activity(
    activity_density: float,
    alpha: float,
    screening_mass: float,
    bridge_conversion: float,
    epsilon: float = 0.0,
) -> float:
    static_source = static_source_amplitude_from_activity(
        activity_density=activity_density,
        bridge_conversion=bridge_conversion,
    )
    mu_value = uniform_screened_mu_from_source(
        static_source=static_source,
        screening_mass=screening_mass,
    )
    return float(tanh_lapse(mu_value, alpha=alpha, epsilon=epsilon))


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


def same_lapse_uniform_wkb_action(
    path_length_m: float,
    lapse: float,
    conserved_energy: float,
    rest_mass: float = 1.0,
) -> float:
    if path_length_m < 0.0:
        raise ValueError("path_length_m must be nonnegative.")
    density = same_lapse_wkb_action_density(
        lapse=np.asarray([lapse], dtype=float),
        conserved_energy=conserved_energy,
        rest_mass=rest_mass,
    )[0]
    return float(path_length_m * density)


def same_lapse_uniform_observables(
    path_length_m: float,
    lapse: float,
    conserved_energy: float,
    rest_mass: float = 1.0,
) -> MassiveModeObservables:
    velocity = same_lapse_group_velocity(
        lapse=np.asarray([lapse], dtype=float),
        conserved_energy=conserved_energy,
        rest_mass=rest_mass,
    )[0]
    turning_lapse = same_lapse_turning_lapse(
        conserved_energy=conserved_energy,
        rest_mass=rest_mass,
    )
    return MassiveModeObservables(
        rest_mass=float(rest_mass),
        conserved_energy=float(conserved_energy),
        turning_lapse=float(turning_lapse),
        minimum_lapse=float(lapse),
        turning_point=0.0 if lapse < turning_lapse else None,
        minimum_group_velocity=float(velocity),
        mean_group_velocity=float(velocity),
        crossing_time_s=same_lapse_uniform_crossing_time_seconds(
            path_length_m=path_length_m,
            lapse=lapse,
            conserved_energy=conserved_energy,
            rest_mass=rest_mass,
        ),
        delay_s=same_lapse_uniform_delay_seconds(
            path_length_m=path_length_m,
            lapse=lapse,
            conserved_energy=conserved_energy,
            rest_mass=rest_mass,
        ),
        wkb_action=same_lapse_uniform_wkb_action(
            path_length_m=path_length_m,
            lapse=lapse,
            conserved_energy=conserved_energy,
            rest_mass=rest_mass,
        ),
    )


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


def same_lapse_species_profile_observables(
    grid_m: np.ndarray,
    lapse: np.ndarray,
    species_parameters: dict[str, tuple[float, float]],
    alpha: float,
    screening_mass: float,
    bridge_conversion: float,
    epsilon: float = 0.0,
) -> dict[str, SpeciesProfileObservables]:
    if bridge_conversion <= 0.0:
        raise ValueError("bridge_conversion must be positive.")
    rows: dict[str, SpeciesProfileObservables] = {}
    for species_name, (rest_mass, conserved_energy) in species_parameters.items():
        observables = same_lapse_profile_observables(
            grid_m=grid_m,
            lapse=lapse,
            conserved_energy=conserved_energy,
            rest_mass=rest_mass,
        )
        required_activity_density = required_activity_density_for_turning(
            conserved_energy=conserved_energy,
            rest_mass=rest_mass,
            alpha=alpha,
            screening_mass=screening_mass,
            bridge_conversion=bridge_conversion,
            epsilon=epsilon,
        )
        rows[str(species_name)] = SpeciesProfileObservables(
            species_name=str(species_name),
            rest_mass=float(rest_mass),
            conserved_energy=float(conserved_energy),
            turning_lapse=float(observables.turning_lapse),
            required_activity_density=float(required_activity_density),
            observables=observables,
        )
    return rows


def theory_native_domain_species_prediction(
    domain: TheoryNativeActivityDomain,
    species: PhysicalSpeciesKinematics,
) -> DomainSpeciesPrediction:
    if domain.activity_density <= 0.0:
        raise ValueError("activity_density must be positive.")
    if domain.path_length_m < 0.0:
        raise ValueError("path_length_m must be nonnegative.")
    bridge_conversion = bridge_conversion_from_theory_native_closure_ratio(
        closure_ratio=domain.closure_ratio,
        screening_length_m=domain.screening_length_m,
        expansion_rate_s_inv=domain.expansion_rate_s_inv,
        tau_planck_multiples=domain.tau_planck_multiples,
    )
    screening_mass = 1.0 / domain.screening_length_m
    uniform_lapse = uniform_lapse_from_activity(
        activity_density=domain.activity_density,
        alpha=domain.alpha,
        screening_mass=screening_mass,
        bridge_conversion=bridge_conversion,
        epsilon=domain.epsilon,
    )
    observables = same_lapse_uniform_observables(
        path_length_m=domain.path_length_m,
        lapse=uniform_lapse,
        conserved_energy=species.conserved_energy_eV,
        rest_mass=species.rest_mass_eV,
    )
    required_activity_density = required_activity_density_for_turning(
        conserved_energy=species.conserved_energy_eV,
        rest_mass=species.rest_mass_eV,
        alpha=domain.alpha,
        screening_mass=screening_mass,
        bridge_conversion=bridge_conversion,
        epsilon=domain.epsilon,
    )
    return DomainSpeciesPrediction(
        domain_name=domain.domain_name,
        species_name=species.species_name,
        bridge_conversion=float(bridge_conversion),
        closure_ratio=float(domain.closure_ratio),
        uniform_lapse=float(uniform_lapse),
        mass_energy_ratio=float(species.mass_energy_ratio()),
        required_activity_density=float(required_activity_density),
        excluded=not math.isfinite(observables.crossing_time_s),
        observables=observables,
    )


def theory_native_domain_species_predictions(
    domain: TheoryNativeActivityDomain,
    species_rows: dict[str, PhysicalSpeciesKinematics],
) -> dict[str, DomainSpeciesPrediction]:
    return {
        str(name): theory_native_domain_species_prediction(domain=domain, species=species)
        for name, species in species_rows.items()
    }


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
