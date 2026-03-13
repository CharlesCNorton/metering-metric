# The Metering Metric

**Duration is not fundamental. It is constituted by measurement.**

A scalar field representing metering density couples to the field Lagrangian as a position-dependent mass term, producing a spacetime geometry in which temporal dynamics are energetically forbidden in regions devoid of metering subsystems. The manifold exists, but nothing elapses.

This is not a restatement of relational time. It is a specific coupling proposal with an abstract theorem development, lattice-model spectral evidence, and falsifiable cosmological predictions.

**Author:** Charles C. Norton

---

## Table of Contents

- [Evidence Boundary](#evidence-boundary)
- [Claim Ladder](#claim-ladder)
- [Framework](#framework)
- [Core Theorems](#core-theorems)
- [The Metering Source](#the-metering-source)
- [Photon Decoupling](#photon-decoupling)
- [Experimental Results](#experimental-results)
  - [Simulator Verification](#simulator-verification)
  - [IBM Quantum Hardware](#ibm-quantum-hardware)
  - [Potential Shape Discrimination](#potential-shape-discrimination)
  - [Coupling Constant Constraints](#coupling-constant-constraints)
- [Metering Field Equation](#metering-field-equation)
- [Cosmological Predictions](#cosmological-predictions)
  - [Void Metering and Horizon Formation](#void-metering-and-horizon-formation)
  - [Neutrino-Photon Timing Delay](#neutrino-photon-timing-delay)
  - [Differential Gravitational Lensing](#differential-gravitational-lensing)
  - [GW170817 Compatibility](#gw170817-compatibility)
- [Observational Proposals](#observational-proposals)
- [Open Questions](#open-questions)
- [Negative Results](#negative-results)
- [Computational Verification](#computational-verification)
- [Status](#status)

---

## Evidence Boundary

This repository mixes several different kinds of support, and they should not be conflated:

- **Coq-formalized:** the root file [`meteringmetric.v`](./meteringmetric.v) currently formalizes the abstract real-analysis skeleton of the framework: proper-time vanishing, tortoise-coordinate divergence, and a WKB-style confinement barrier theorem under stated hypotheses on the activation function and metering profile.
- **Computationally checked:** the lattice, PDE, and cosmology sections are supported by numerical experiments and symbolic derivations, not by formal proof.
- **Hardware-backed:** the IBM Quantum section is evidence for discrete bound-state structure in the tested lattice model, not a direct measurement of the physical spacetime coupling constant.
- **Measurement-anchored laboratory benchmark:** [`laboratory_bounds.py`](./laboratory_bounds.py) now maps the repository's current Casimir-scale anomaly benchmark onto published precision Casimir-pressure nulls, a first active-vs-inert differential design target, a first geometry-dependent planar slab signal model sourced by the screened metering field equation, and a first explicit source-bridge report. Under the present benchmark model, the strongest executed laboratory null in the repository is a Casimir bound of $\alpha \leq 6.79 \times 10^{-3}$, with canonical reports in [`laboratory_reports/casimir-benchmark-report.json`](./laboratory_reports/casimir-benchmark-report.json), [`laboratory_reports/casimir-differential-design.json`](./laboratory_reports/casimir-differential-design.json), [`laboratory_reports/casimir-planar-slab-geometry.json`](./laboratory_reports/casimir-planar-slab-geometry.json), and [`laboratory_reports/casimir-source-bridge-report.json`](./laboratory_reports/casimir-source-bridge-report.json). In the current slab model, the interaction law is exact, the naive bridge is not viable, and the first viable bridge family is now a Planck-tick occupancy reduction from dynamic activity density to static screened source.
- **Archival observational program:** the cluster-lensing section currently runs on two HFF fields, a hardened radial-median-bandpass residual, BUFFALO v2 photometric catalogs, public MUSE spectroscopic catalogs, deeper external spectroscopy for membership hardening, public `MUSE-DEEP` core-cube extractions, and 22 independent lens-model families. The sign-flip morphology survives those checks. The most stable explanatory axes are local convergence and stellar mass. Generic spectroscopic structure, line-ratio proxies, first strong-line metallicity estimators, and the current public matched-depth core-cube route do not yield a same-sign replicated abundance-style proxy across both clusters under the present controls.
- **Open empirical question:** the cosmological interpretation, observer-detection ideas, and any claim that the measured morphology tracks decoherence-rate proxies beyond mass and local convergence remain open. The observational proxy lane is exhausted on current public data. Two later observational pursuits remain well-defined: redder symmetric spectroscopy for both clusters with usable `Halpha`, `[N II]`, `Hbeta`, and `[O III]` coverage, and direct per-galaxy abundance tables with one-to-one mapping into the HFF member catalogs. Until such data exist, the active focus shifts to deriving the 3+1D lensing observable from the field equations.

---

## Claim Ladder

The repository now supports a cleaner claim hierarchy:

- **Established within the current project:** the mass-gap mechanism is the working core. The abstract theorem spine in [`meteringmetric.v`](./meteringmetric.v) supports horizon formation, confinement-barrier structure, and proper-time vanishing in the present theorem-level scope. The lattice and hardware program supports discrete bound-state structure for the tested $1/N^2$ potential. The archival lensing program supports a reproducible sign-flip-like residual morphology, with local convergence and stellar mass as the stable explanatory axes on current public data.
- **Supported but not finished:** a viable 3+1D Einstein-backreaction branch exists in [`three_plus_one_lensing.py`](./three_plus_one_lensing.py). In that branch, raw convergence remains nonnegative in the tested cluster-style maps, but sign-changing residual morphology appears in multi-component, HFF-member, and calibrated member-plus-envelope geometries. Direct archive comparison improves under that calibration, and a rigid-alignment audit shows that part of the remaining mismatch is registration-sensitive.
- **Open:** a full 3+1D observable derivation is still needed. The current theory program has not yet turned the viable branch into a quantitatively competitive field-level match to the archival residual maps. A continuum finite-dimensionality theorem remains open. A direct first-principles physical measurement of the spacetime coupling $\alpha$ remains open, even though the repository now contains a measurement-anchored Casimir benchmark null bound, a differential design target, an exact planar screened-slab interaction law, and a first viable source-bridge family. On the laboratory side, the remaining gap is no longer force scaling and no longer an unconstrained normalization constant; it is the covariant fixation of the occupancy bridge parameters $(\eta_J, \tau_p, L_\perp)$. On the lensing side, the present Einstein scaffold fixes only a composite backreaction scale rather than the source bridge parameters separately. On the observational side, the next clean empirical step would require redder symmetric spectroscopy or direct abundance tables for both clusters.
- **Not supported yet:** the repository does **not** currently support the stronger claim that the archival cluster morphology tracks decoherence-rate or abundance-style proxies beyond mass and local convergence. It does **not** yet support a raw 3+1D negative-$\kappa$ prediction in the viable branch. It does **not** yet support treating the explicit photon-sector branch as the mainline completion.

This ladder is the current intended reading of the project: the core metering mechanism remains alive, the observational proxy claim has not cleared the bar, and the active bridge to reality is the 3+1D Einstein-backreaction lensing program.

---

## Framework

### Intellectual Context

The framework draws on five established programs, each of which stops short of the specific mechanism formalized here:

- **Wheeler-DeWitt (1967):** The wavefunction of the universe satisfies H|Psi> = 0. No time parameter appears. Time is absent at the fundamental level.
- **Page-Wootters (1983):** A timeless universe contains time-evolving subsystems when conditioned on a clock degree of freedom. Time is created by conditioning, not discovered.
- **Rovelli (1996, 2018):** Duration is a relation between systems. The thermal time hypothesis derives time flow from thermodynamic structure.
- **Barbour (1999):** Time can be eliminated from the fundamental description, replaced by relations between configurations.
- **Bridgman (1927):** A quantity without a measurement procedure is undefined. Unmeasured duration is not duration.

None of these produces the specific mass-term coupling proposed here. This work develops that proposal into a concrete analytic and computational program.

### The Metric

In 1+1 dimensions:

$$ds^2 = -N(\mu(x))^2 \, dt^2 + dx^2$$

where:
- $\mu(x)$ is the **metering density**: the local density of systems capable of recording state transitions
- $N(\mu) = \tanh(\alpha \cdot \mu)$ is the **lapse function**, monotone increasing with $N(0) = 0$ and $N(\infty) \to 1$
- $\alpha$ is the coupling constant

The coupling enters through the **mass term**. For a massive field $\phi$:

$$\mathcal{L}_{\text{massive}} = \frac{1}{2N^2}(\partial_t \phi)^2 - \frac{1}{2}(\partial_x \phi)^2 - \frac{1}{2N^2}\phi^2$$

The effective mass is $m_{\text{eff}}(x) = 1/N(x)$, which diverges as $\mu \to 0$. Temporal oscillations become energetically forbidden. The field does not "slow down" in unmetered regions -- its temporal degree of freedom freezes.

### Why Mass, Not Geometry

The original approach (v1) attempted to implement conditional duration via the lapse function in the geodesic equation. This failed: geodesic equations enforce their own normalization, so worldlines accumulate proper time regardless of $N$. The mass-gap formulation succeeds because it makes temporal oscillation infinitely costly, not geometrically absent. The mechanism is energetic, not kinematic.

---

## Core Theorems

Three structural results motivate the framework for **any** monotone activation function $f$ (tanh, sigmoid, arctan, erf) and **any** metering profile $\mu$ that decays to zero (Gaussian, top-hat, power-law, double-Gaussian). Their abstract analytic skeleton is formalized in Coq, while the dimensional extensions discussed here are computational rather than machine-checked.

### Theorem 1: Horizon Formation

The tortoise coordinate integral

$$x^*(L) = \int_0^L \frac{dx}{N(x)}$$

diverges as $L \to \infty$. The metering boundary is a true horizon: signals from the metered region take infinite coordinate time to reach the unmetered region.

**Proof sketch:** Since $\mu$ decays to zero, $N(\mu(x)) \to 0$ for large $|x|$. The integrand $1/N$ is unbounded, and the integral diverges.

This has been verified computationally across 4 activation functions, 4 metering profiles, and 3 spatial dimensions.

### Theorem 2: Temporal Mode Confinement

The effective potential $V_{\text{eff}}(x) = 1/N(x)^2$ confines temporal oscillations:

- $V_{\text{eff}} \to \infty$ as $\mu \to 0$ (the potential is confining)
- The abstract theorem spine proves a WKB-style barrier under the stated decay and positivity hypotheses
- In the representative lattice and continuum truncations studied here, the temporal spectrum is discrete and bound-state dominated
- In the representative WKB estimate used here, the tunneling probability is numerically negligible: $\exp(-3218)$
- The representative numerical mode counts are: ~69 radial modes (1+1D), ~600 modes (1000-site continuum), ~17 bits (3+1D with angular modes)

These finite counts are numerical results for the tested truncations and dimensional models. A continuum finite-dimensionality theorem remains open.

The current same-lapse massive-sector workbench now gives explicit turning thresholds, barrier actions, and delay integrals rather than only qualitative confinement statements.

### Theorem 3: Proper Time Vanishing

For a static observer at position $x$:

$$\tau(x, T) = N(x) \cdot T \to 0 \quad \text{as} \quad \mu(x) \to 0$$

In the current formal development, if a prescribed path $\gamma(t)$ lies entirely in a meterless region on an interval $[t_1, t_2]$, then the lapse integral over that interval is exactly zero.

This is the present theorem-level scope. Extending it to a full dynamical statement for general worldlines is a separate derivation.

### Formal Verification

A Coq development is now present in the root file [`meteringmetric.v`](./meteringmetric.v). It compiles locally against Coquelicot and establishes the abstract theorem spine used in this repository:

- proper-time vanishing in meterless regions
- unbounded tortoise-coordinate growth
- a WKB-style confinement barrier result under the stated decay and positivity hypotheses

What it does **not** yet provide is a full formalization of the README's physical superstructure. The numerical dimension-specific counts, PDE dynamics, hardware interpretation, and cosmological phenomenology remain computational or interpretive results rather than machine-checked theorems.

---

## The Metering Source

The metering field equation has a source term $J(x)$ that determines where $\mu$ is nonzero. Three candidate criteria were evaluated:

| Criterion | $J$ in void | $J$ in diffuse gas | $J$ in lab | Horizons? |
|-----------|-------------|---------------------|------------|-----------|
| Entropy production rate | ~0 | small but > 0 | large | Marginal |
| Structured record capacity | 0 | ~0 | large | Yes |
| Decoherence rate density | 0 | tiny | large | Yes |

### Recommended: Decoherence Activity Density

Define the dynamic metering-activity density

$$R(x) = \sum_{\text{massive } i} n_i(x) \cdot \gamma_{D,i}(x)$$

where $n_i$ is the number density of massive species $i$ and $\gamma_{D,i}$ is the local decoherence rate. The current covariant mainline does not identify the screened source with raw $R$ directly. It introduces an occupancy density $Q$ carried by the preferred frame $u^a$ and transported by

$$\nabla_a(Q\,u^a) = R - \frac{Q}{\tau_p}.$$

The screened source entering the static $\mu$ equation is then

$$J(x) = \eta_J \frac{Q(x)}{L_\perp(x)^2}.$$

In a quasistatic local rest frame with expansion scalar $\theta = \nabla_a u^a$ and negligible convective derivative $u^a\nabla_a Q$, this reduces to

$$
Q = \frac{\tau_p}{1+\theta \tau_p}R,
\qquad
J(x) = \kappa_J^{\mathrm{eff}}(x)\,R(x),
\qquad
\kappa_J^{\mathrm{eff}}(x) = \eta_J \frac{\tau_p(x)}{L_\perp(x)^2\left(1+\theta(x)\tau_p(x)\right)}.
$$

Here $\tau_p$ is the persistence time of a metering event in the static source and $L_\perp$ is the transverse coarse-graining scale of that reduction.

**Why this criterion:**

1. **The dynamic source is fixed by local activity.** The entropy production criterion requires an additional coupling $\beta$; the structured record criterion requires an arbitrary persistence threshold $\tau_{\text{threshold}}$. In the decoherence-activity picture, the dynamical source density itself is fixed by local rate density, and the remaining bridge to the static screened source is now explicit and testable.
2. **Consistent with photon decoupling.** Photons don't decohere in vacuum. The CMB photon bath is already at thermal equilibrium. Only massive degrees of freedom under environmental monitoring contribute.
3. **Physically grounded.** Decoherence is measurement in the decoherence program. The framework claims duration requires measurement. The source of metering should be the rate of measurement events.
4. **Correct hierarchy:**

| Environment | $\gamma_D$ (s$^{-1}$) | $\rho_{\text{states}}$ (m$^{-3}$) | $R = \sum n_i\gamma_{D,i}$ (m$^{-3}$ s$^{-1}$) |
|-------------|------------------------|--------------------------------------|---------------------------|
| Lab (condensed matter) | $10^{12}$ | $10^{28}$ | $10^{40}$ |
| Stellar core | $10^{15}$ | $10^{31}$ | $10^{46}$ |
| Molecular cloud | $10^{3}$ | $10^{10}$ | $10^{13}$ |
| Diffuse ISM | $10^{-5}$ | $10^{6}$ | $10^{1}$ |
| IGM | $10^{-15}$ | $1$ | $10^{-15}$ |
| Cosmic void (baryons) | $10^{-20}$ | $0.1$ | $10^{-21}$ |
| CMB photon field | $0$ | $0$ (massless) | $0$ |

The laboratory workbench now fixes the bridge scale at the reference slab benchmark point. In [`laboratory_reports/casimir-source-bridge-report.json`](./laboratory_reports/casimir-source-bridge-report.json), the exact slab law selects

$$\kappa_J \approx 1.90 \times 10^{-32}\ \text{s m}^{-2}.$$

That means the naive raw-rate identification is out. The first viable bridge family is instead a Planck-tick occupancy reduction,

$$J(x) = \eta_J \frac{t_P}{L_\perp(x)^2} R(x).$$

For the current slab geometry, taking $L_\perp = \ell_s = 0.2\,\mu\text{m}$ requires

$$\eta_J \approx 1.41 \times 10^{-2},$$

while taking $\eta_J = 1$ instead selects

$$L_\perp \approx 1.68\,\mu\text{m}.$$

So the bridge is no longer an arbitrary missing constant. It has been reduced to a narrow occupancy-style family that now needs to be fixed by the covariant theory.

The natural closure of that family uses the two scales already present in the active theory:

$$
\tau_p = \chi_\tau\, t_P,
\qquad
L_\perp = \chi_L\,\ell_s,
\qquad
\ell_s = m_\mu^{-1}
$$

in units with $c=1$, or $\ell_s$ equal to the physical screening length in SI units. Then

$$
\kappa_J^{\mathrm{eff}}
=
\eta_J\,
\frac{\chi_\tau t_P}{\chi_L^2 \ell_s^2 \left(1+\theta \chi_\tau t_P\right)}.
$$

So once the source reduction is expressed in the natural metering scales, the bridge freedom collapses from three dimensional quantities to one dimensionless closure ratio

$$
\Xi_J = \eta_J \frac{\chi_\tau}{\chi_L^2}.
$$

At the laboratory benchmark point with $\ell_s = 0.2\,\mu\mathrm{m}$ and negligible expansion correction, the exact slab law fixes

$$
\Xi_J \approx 1.41 \times 10^{-2}.
$$

In the minimal natural reduction, that ratio has a direct theory meaning rather than a purely phenomenological one. If a metering event occupies one transverse screening cell of area $\ell_s^2$ for $\chi_\tau t_P$ and the source reduction counts only the occupied fraction of those cells, then

$$
\Xi_J
=
\eta_J \frac{\chi_\tau}{\chi_L^2}
$$

is exactly the dimensionless screening-cell occupancy fraction carried into the static screened source. The current benchmark value

$$
\Xi_J \approx 1.41 \times 10^{-2}
$$

therefore reads as an occupancy fraction of about `1.4%` per natural screening cell in the minimal closure, not as an arbitrary dimensional fudge factor.

That is the current strict mainline reading of the bridge problem. The theory no longer needs to determine three unrelated dimensional inputs. It needs to determine the dimensionless closure ratio $\Xi_J$ and the circumstances under which the natural-scale identification $(\tau_p,L_\perp) = (\chi_\tau t_P,\chi_L \ell_s)$ is correct.

The new cross-domain audit in [`laboratory_reports/source-closure-crossdomain.json`](./laboratory_reports/source-closure-crossdomain.json) now shows exactly where the literal one-cell reading breaks. If one takes the Planck-tick screening-cell picture at face value and computes the naive occupied-cell fraction

$$
f_{\mathrm{cell}} \sim R\,t_P\,\ell_s^3,
$$

then at the laboratory benchmark point `R_{\mathrm{lab}} \sim 10^{40}\ \mathrm{m^{-3}s^{-1}}` with `\ell_s = 0.2\,\mu\mathrm{m}` one gets only

$$
f_{\mathrm{cell,lab}} \approx 4.31 \times 10^{-24},
$$

so matching the current benchmark closure ratio would require a collective enhancement

$$
\Xi_J / f_{\mathrm{cell,lab}} \approx 3.27 \times 10^{21}.
$$

At the current lensing reference source amplitude `10^{-3}` on that same bridge, the implied activity scale `R \approx 5.26 \times 10^{28}\ \mathrm{m^{-3}s^{-1}}` gives

$$
f_{\mathrm{cell,lensing}} \approx 2.27 \times 10^{-35},
$$

so the required enhancement rises to

$$
\Xi_J / f_{\mathrm{cell,lensing}} \approx 6.21 \times 10^{32}.
$$

But in the `1 Mpc` reference-void branch with `R_{\mathrm{void}} \sim 10^{-21}\ \mathrm{m^{-3}s^{-1}}`, the same literal rule oversaturates:

$$
f_{\mathrm{cell,void}} \approx 1.58 \times 10^3,
\qquad
\Xi_J / f_{\mathrm{cell,void}} \approx 8.90 \times 10^{-6}.
$$

So the literal Planck-tick single-screening-cell picture is not a universal closure of the mainline. What survives is the covariant transport law together with a renormalized effective reduction represented by $\Xi_J$. The open bridge problem is now narrower and sharper: derive that renormalized reduction from the source sector itself, rather than treating it as a benchmark-matched constant.

The same bridge law now ties the three active domains together:

$$
J(x) = \kappa_J^{\mathrm{eff}}(x)\,R(x), \qquad
S_{\mathrm{eff}}(x) \propto \texttt{gravity\_scale}\,\texttt{density\_scale}\,\left(\kappa_J^{\mathrm{eff}}(x)\right)^2 R(x)^2.
$$

So once $(\eta_J,\tau_p,L_\perp)$ are fixed, the domain ratios are no longer free:

$$
\frac{J_A}{J_B} = \frac{\kappa_{J,A}}{\kappa_{J,B}} \frac{R_A}{R_B},
$$

and in the shared-bridge specialization $\kappa_{J,A} = \kappa_{J,B}$ this reduces to $J_A/J_B = R_A/R_B$. The laboratory program fixes the benchmark size of $\kappa_J$, the void calculations use the same law to determine which activity contrasts can sustain horizons, and the lensing branch uses the same law to set the source amplitude entering the screened Einstein backreaction.

That same bridge now gives a concrete scale for the present lensing scaffold. If the current reference Einstein calibration uses source amplitude `10^{-3}` and the laboratory bridge value above, then the implied activity density is

$$
R_{\mathrm{lens,ref}} = \frac{10^{-3}}{\kappa_J} \approx 5.26 \times 10^{28}\ \mathrm{m}^{-3}\mathrm{s}^{-1}.
$$

That number is not yet a full physical extraction of cluster metering activity, because the current lensing scaffold still determines only the composite backreaction scale $\texttt{gravity\_scale}\,\texttt{density\_scale}\,\texttt{source\_amplitude}^2$. But it does fix the bridge-side scale that any completed cross-domain theory has to respect.

For a Planck-tick bridge, the explicit expansion correction is tiny in the present universe. Using $\theta \sim 3H_0$ with $H_0 \approx 2.2 \times 10^{-18}\ \mathrm{s}^{-1}$ gives

$$\theta t_P \approx 3.56 \times 10^{-61},$$

so the lab- and cluster-scale bridge is effectively the $\theta = 0$ limit unless $\tau_p$ is vastly larger than $t_P$.

5. **Reduces to entropy production in the thermodynamic limit** via the fluctuation-dissipation theorem.

---

## Photon Decoupling

### The Problem

If all fields couple to the metering metric, then photon propagation speed is $N(x)$, and light cannot traverse cosmic voids ($N \to 0$). We observe photons from 13 billion light-years away through voids.

### The Resolution

The coupling enters through the **mass term**. Photons are massless -- there is no mass term to diverge. They propagate on the bare metric $ds^2 = -dt^2 + dx^2$ at $c$ everywhere.

The coupling hierarchy:

1. **The metering field $\mu$ itself** couples to $N(\mu)$: propagation speed $c_\mu = N$. Self-reinforcing horizons.
2. **Massive matter fields** couple to $N(\mu)$: effective mass $m_{\text{eff}} = 1/N$ diverges in voids. Temporal dynamics freeze.
3. **Massless radiation** (photons, gravitons, gluons) couples to the bare metric: propagates at $c$ regardless of $\mu$.

This is a species-dependent metric -- a bimetric structure. The conditional temporal metric already breaks the equivalence principle by coupling the metric to metering density rather than stress-energy. Photon decoupling is a refinement of an existing violation, not a new one.

### What This Preserves

Temporal horizons, mode confinement, bounded temporal information, self-reinforcing horizons, geometric SETI signatures, and WEC violation structure are all unchanged within the model. Photon propagation through voids is restored. Under this mass-term-only coupling, the GW170817 constraint is satisfied.

### What This Costs

Equivalence principle violation for massive vs. massless fields. Bimetric and disformal coupling theories are established in modified gravity (Bekenstein's TeVeS, massive gravity). The mathematical framework exists.

---

## Experimental Results

### Simulator Verification

#### 6-Qubit Density Matrix (Unity laptop)

6-site lattice, Gaussian metering profile ($\sigma = 1.5$), dephasing Kraus channels modeling spatially varying measurement strength. Single-excitation subspace.

Predicted bound-state energies:

| $E_0$ | $E_1$ | $E_2$ | $E_3$ | $E_4$ | $E_5$ |
|-------|-------|-------|-------|-------|-------|
| 0.487 | 2.057 | 3.683 | 4.113 | 16.829 | 16.829 |

Survival probability shows oscillatory revival structure -- **not** monotonic Zeno decay. This is the signature of a discrete bound-state spectrum.

Fourier spectral peaks:

| Extracted $\omega$ | Predicted transition | Relative error |
|--------------------|----------------------|----------------|
| 1.745 | $E_2 - E_1 = 1.626$ | 7.3% |
| 3.491 | $E_3 - E_0 = 3.626$ | 3.7% |

Error limited by 48-point time series resolution ($\Delta\omega = 0.87$).

#### 1000-Site Continuum Limit (SAURON, i9-13900KF + RTX 6000 Ada)

1000 x 1000 density matrix, single-excitation subspace. Gaussian metering ($\sigma = 100$). 400 time points, 8 Trotter steps each. Eigendecomposition in 0.31s, full evolution in 1316s.

- **600 bound states** below the potential barrier
- Spacing ratio $(E_2 - E_1)/(E_1 - E_0) = 1.0069$ (nearly harmonic at the bottom, expected for smooth confining potential)
- **Five Fourier peaks, all matching predicted eigenvalue differences:**

| Extracted $\omega$ | Predicted transition | Relative error |
|--------------------|----------------------|----------------|
| 0.7854 | $E_{169} - E_{129} = 0.7854$ | 0.0016% |
| 1.5708 | $E_{131} - E_{66} = 1.5699$ | 0.058% |
| 2.3562 | $E_{145} - E_{46} = 2.3573$ | 0.047% |
| 3.1416 | $E_{167} - E_{30} = 3.1406$ | 0.033% |
| 3.9270 | $E_{177} - E_0 = 3.9269$ | 0.003% |

Every peak matches to better than **0.06%**. This is two orders of magnitude tighter than the 6-qubit result and is consistent with convergence toward a continuum bound-state structure for the tested $1/N^2$ potential.

The continuum statement is now narrower and stronger than the older finite-count wording. For the one-dimensional Dirichlet half-line operator

$$
H_i = -\partial_x^2 + V^{(i)}_{\mathrm{eff}}(x), \qquad
V^{(i)}_{\mathrm{eff}}(x) = \frac{m_i^2}{N(x)^2},
$$

the repository's formal spine now proves the key structural hypotheses on the active branch: $V^{(i)}_{\mathrm{eff}}$ is positive and diverges as $x \to +\infty$. The intended continuum claim is therefore not a finite-dimensional Hilbert space claim. It is a discrete-spectrum claim for the confining massive-sector operator.

The proof on the half-line is standard but can be stated directly. Let

$$
H_i = -\partial_x^2 + V^{(i)}_{\mathrm{eff}}(x)
$$

on $L^2([0,\infty))$ with Dirichlet boundary at $x=0$, and assume $V^{(i)}_{\mathrm{eff}}$ is continuous, nonnegative, and diverges to $+\infty$. Define the quadratic-form domain

$$
\mathcal D(q_i)
=
\left\{
u \in H_0^1([0,\infty)) :
\int_0^\infty V^{(i)}_{\mathrm{eff}}(x)\,|u(x)|^2\,dx < \infty
\right\}
$$

with

$$
q_i[u] = \int_0^\infty \left(|u'(x)|^2 + V^{(i)}_{\mathrm{eff}}(x)|u(x)|^2\right)\,dx.
$$

Fix any sequence bounded in the form norm. Because $V^{(i)}_{\mathrm{eff}}(x)\to+\infty$, for every $M>0$ there exists $L_M$ such that $V^{(i)}_{\mathrm{eff}}(x)\ge M$ on $[L_M,\infty)$. Therefore

$$
\int_{L_M}^\infty |u(x)|^2\,dx
\le
\frac{1}{M}\int_{L_M}^\infty V^{(i)}_{\mathrm{eff}}(x)|u(x)|^2\,dx
\le
\frac{1}{M} q_i[u].
$$

So a form-bounded sequence has uniformly small $L^2$ tail once $M$ is large. On the finite interval $[0,L_M]$, the same sequence is bounded in $H_0^1([0,L_M])$, and the one-dimensional Rellich theorem gives precompactness in $L^2([0,L_M])$. Combining the compact interior piece with the uniformly small tail shows that the embedding

$$
\mathcal D(q_i) \hookrightarrow L^2([0,\infty))
$$

is compact. Hence the Friedrichs extension of $H_i$ has compact resolvent, and therefore a purely discrete spectrum of finite-multiplicity eigenvalues accumulating only at $+\infty$.

That is the continuum result the project uses. The current finite-difference half-line workbench now matches that reading numerically as well: for the reference confining branch with $N(x)=\tanh(e^{-x})$, the first five Dirichlet eigenvalues stabilize to displayed precision under box enlargement from length `6` to `10`, namely `9.4615`, `21.2414`, `35.2478`, `51.1856`, and `68.8673`.

### IBM Quantum Hardware

Three runs on `ibm_torino` (Qiskit 2.3.0, SamplerV2), 6-qubit lattice, 4000 shots per circuit.

**Run 1** (14 circuits, $t = 0.15$--$2.10$, Job `d6e3gdp54hss73badhbg`):
Survival probability descends from 0.87 to 0.13, then revives to 0.33. Dominant Fourier peak $\omega = 2.99$ matches predicted $E_2 - E_0 = 3.20$ at 6.4% error.

**Run 2** (20 circuits, $t = 0.05$--$1.00$, Job `d6e3ih154hss73badjt0`):
Dense sampling of the descent phase. Hardware tracks noiseless simulator within 3--8% for circuits with $\leq 50$ two-qubit gates.

**Run 3** (20 circuits, $t = 0.30$--$2.50$, Job `d6e3l9vg4t5c7387b93g`):
Independent confirmation. Minimum $P = 0.15$ at $t = 1.23$, revival to $P = 0.28$ at $t = 1.57$, second descent to $P = 0.16$ at $t = 2.38$.

The oscillatory revival structure -- descent, minimum, revival, second descent -- is **reproduced across two independent submissions** with different circuit constructions. Within the tested lattice model, this argues against generic monotonic Zeno freezing and supports a bound-state spectral interpretation.

### Potential Shape Discrimination

Chi-squared fitting of hardware data against five potential shapes, each matched in depth and width:

| Potential | Run 1 $\chi^2$ | Run 3 $\chi^2$ | Rank |
|-----------|-----------------|-----------------|------|
| $1/N^2$ (metering) | 1113 | 703 | **1st** |
| Quartic ($x^4$) | 1991 | 1481 | 2nd |
| Harmonic ($x^2$) | 6902 | 10344 | 3rd |
| Square well | 7277 | 16709 | 5th |
| Linear ($|x|$) | 7624 | 14769 | 4th |

**Discrimination ratios:**
- $1/N^2$ vs. quartic: **1.8--2.1x** (nearest competitor by $\chi^2$ in this test)
- $1/N^2$ vs. harmonic: **6--15x**

In this 6-site hardware comparison, the $1/N^2$ potential gives the best fit among the tested alternatives. The key discriminant is the potential's behavior at intermediate sites: $1/N^2$ gives $V = 3.41$ where quartic gives $V = 3.77$ and harmonic gives $V = 7.21$. The steep walls of the $1/N^2$ potential (rapid transition from low to high $V$) distinguish it from smoother shapes.

### Coupling Constant Constraints

Four current handles on the coupling constant $\alpha$:

**1. Hardware spectral fit:** $\alpha = 1.165$ minimizes $\chi^2$ against `ibm_torino` data (combined Runs 1 + 3).

**2. Internal weak-energy benchmark:** The metric requires negative energy density inside the metered region at ~16x the Casimir scale when $\alpha = 1$. Requiring $|T_{00}| < E_{\text{Casimir}}$ gives $\alpha < 0.133$.

**3. Executed Casimir benchmark null:** The repository now carries an explicit laboratory workbench in [`laboratory_bounds.py`](./laboratory_bounds.py), with the canonical null report in [`laboratory_reports/casimir-benchmark-report.json`](./laboratory_reports/casimir-benchmark-report.json). Using the same Casimir-scale benchmark amplitude as above, together with quadratic weak-coupling scaling, the strongest current pressure-null anchor in the report is the Decca et al. 2003 parallel-plate measurement at `0.2 um`, which gives

$$\alpha \leq 6.79 \times 10^{-3}.$$

The secondary Decca et al. 2007 relative-error anchor gives

$$\alpha \leq 1.09 \times 10^{-2}.$$

At the `0.2 um` reference separation used in the report, this benchmark model implies:

- $\alpha = 10^{-2}$ would require an anomalous pressure of `1.30 mPa`
- $\alpha = 10^{-3}$ would require an anomalous pressure of `13.0 uPa`

So the first hard laboratory cutoff now exists in the repository, and it is already stronger than the older internal WEC consistency cutoff.

The same workbench now also carries the first active-vs-inert differential design report in [`laboratory_reports/casimir-differential-design.json`](./laboratory_reports/casimir-differential-design.json). At the current reference design point:

- separation `0.2 um`
- active area `1 mm^2`
- target `\alpha = 10^{-3}`

the benchmark differential is `13.0 uPa` or `13.0 pN`, with required common-mode rejection `6.25 x 10^4`. If uncompensated gap drift were allowed to fake the whole signal, the corresponding gap-stability scale would be `0.8 pm`.

**4. Dimensional analysis:** In SI units, $\alpha_{\text{SI}} \sim \ell_P^3 \sim 4 \times 10^{-105}$ m$^3$/bit. In a laboratory ($\mu \sim 10^{26}$ bits/m$^3$): $\alpha_{\text{SI}} \cdot \mu \sim 4 \times 10^{-79}$.

**Interpretation:** The hardware test supports the lattice-model claim that a $1/N^2$ potential can produce the observed bound-state spectral structure. The Casimir workbench now adds the first executed measurement-anchored laboratory null bound, the first concrete differential design target, the first geometry-dependent slab-field signal model in the repository, and the first explicit source-bridge report. The planar interaction shape is now exact within the linear screened slab model, and the finite-domain solver agrees with that closed form at the `10^{-5}` level. The raw dynamic activity density $R$ is not itself the static screened source, and the naive bridge is already excluded at the current benchmark scale. The current live bridge family is instead occupancy-style, with an occupancy density $Q$ obeying $\nabla_a(Q u^a) = R - Q/\tau_p$ and static source $J = \eta_J Q/L_\perp^2$; in the quasistatic limit that reduces to $J = \eta_J \tau_p R / (L_\perp^2(1+\theta \tau_p))$. A direct physical extraction of spacetime $\alpha$ therefore remains open, but it is now constrained to a much narrower source-bridge problem. The active separation of lattice-fit, continuum, and laboratory couplings is now built directly into the repository reading itself: $\alpha_{\text{lattice}} = 1.165$, the internal continuum consistency cutoff is `0.133`, and the strongest current benchmark laboratory null is `6.79 x 10^{-3}`.

---

## Metering Field Equation

### Action and Dynamics

The active completion is now a preferred-frame covariant theory with four ingredients:

- a spacetime metric $g_{ab}$
- a metering scalar $\mu$
- a unit timelike frame $u^a$
- an environment/source sector $\Psi$ that coarse-grains to $R$
- massive matter fields $\Phi_i$ plus a photon field $A_a$

Write

$$N(\mu) = \epsilon + (1-\epsilon)\tanh(\alpha \mu), \qquad Z(\mu) = N(\mu)^{-2}.$$

Then the mainline action is

$$
S = \int d^4x \sqrt{-g}\left[
\frac{R[g]}{16\pi G}
- \frac12 \nabla_a\mu \nabla^a\mu
- \frac12 m_\mu^2 \mu^2
+ \lambda (u_a u^a + 1)
\right]
+ S_{\mathrm{src}}[\Psi,\mu; g,u]
+ \sum_i S_i[\Phi_i; g,\mu]
+ S_\gamma[A;g],
$$

with

$$
S_i = -\frac12 \int d^4x \sqrt{-g}\left[
 g^{ab}\nabla_a\Phi_i\nabla_b\Phi_i
+ Z(\mu)\,m_i^2 \Phi_i^2
\right],
$$

and

$$
S_\gamma[A;g] = -\frac14 \int d^4x \sqrt{-g}\,F_{ab}F^{ab}.
$$

The minimal closed source sector used by the current mainline introduces an occupancy scalar $Q$ together with a multiplier field $\Lambda$ enforcing the source-transport law:

$$
S_{\mathrm{src}}[\Psi,\mu;g,u]
= \int d^4x \sqrt{-g}\left[
\mathcal{L}_{\mathrm{occ}}(\Psi,\nabla\Psi;g,u)
+ \Lambda\left(\nabla_a(Q u^a) - R(\Psi,u,g) + \frac{Q}{\tau_p(\Psi,u,g)}\right)
+ \mu\,\eta_J \frac{Q}{L_\perp(\Psi,u,g)^2}
\right].
$$

This is the strict mainline reading of the repository:

- the massive sector couples through the mass term, with effective potential $m_i^2/N(\mu)^2$
- photons propagate on the metric $g_{ab}$ itself
- lensing and redshift therefore come from the same metric once the Einstein equations are solved

The source is therefore the reduced static screened scalar built directly into the source/environment sector,

$$
J = \eta_J \frac{Q}{L_\perp^2}, \qquad
\nabla_a(Q u^a) = R - \frac{Q}{\tau_p}, \qquad
R = \sum_i n_i\,\gamma_{D,i},
$$

with $R$, $\tau_p$, and $L_\perp$ treated as source-sector scalars in the local preferred rest frame. In the quasistatic reduction this closes to

$$
Q = \frac{\tau_p}{1+\theta\tau_p}R, \qquad
J = \eta_J \frac{\tau_p}{L_\perp^2(1+\theta\tau_p)}R,
$$

where $\theta = \nabla_a u^a$.

The current mainline closure specializes those source-sector scalars to the natural metering scales,

$$
\tau_p = \chi_\tau t_P,
\qquad
L_\perp = \chi_L \ell_s,
\qquad
\ell_s = m_\mu^{-1}
$$

in units with $c=1$, so the static screened source becomes

$$
J
=
\eta_J
\frac{\chi_\tau t_P}{\chi_L^2 \ell_s^2 (1+\theta \chi_\tau t_P)}\,R
=
\kappa_J^{\mathrm{eff}} R.
$$

That is the covariant closure now carried by the repository: the theory fixes a dimensionless source-reduction ratio rather than three unrelated dimensional bridge parameters.

Within the current static Einstein workbench, one further scaling fact is exact. For fixed source-shape parameters, the screened field equation is linear in the source amplitude, so $\mu \propto \texttt{source\_amplitude}$. The scalar stress tensor is quadratic in $\mu$ and $\partial\mu$, so

$$
T^{(\mu)}_{ab} \propto \texttt{density\_scale}\,\texttt{source\_amplitude}^2.
$$

The backreacted photon metric therefore depends on the current scaffold only through the composite scale

$$
S_{\mathrm{eff}} = \texttt{gravity\_scale}\,\texttt{density\_scale}\,\texttt{source\_amplitude}^2
$$

once the geometric shape parameters are held fixed. That exact scaling is now confirmed numerically in the HFF comparison workbench as well.

The corresponding field equations are:

$$
G_{ab} = 8\pi G\left(T^{(\mu)}_{ab} + T^{(u)}_{ab} + T^{(\mathrm{src})}_{ab} + \sum_i T^{(i)}_{ab} + T^{(\gamma)}_{ab}\right),
$$

$$
\Box\mu - m_\mu^2 \mu + J - \frac12 Z'(\mu)\sum_i m_i^2 \Phi_i^2 = 0,
$$

$$
\nabla_b F^{ab} = 0,
$$

$$
\nabla_a(Q u^a) = R - \frac{Q}{\tau_p},
$$

together with the unit constraint $u_a u^a = -1$, the source-sector Euler-Lagrange equations from $S_{\mathrm{src}}$, and the matter Euler-Lagrange equations from $S_i$. In particular, varying the multiplier field gives the transport law itself, while varying $Q$ closes the source sector through

$$
u^a \nabla_a \Lambda - \frac{\Lambda}{\tau_p} + \eta_J \frac{\mu}{L_\perp^2} + \frac{\partial \mathcal{L}_{\mathrm{occ}}}{\partial Q} = 0.
$$

Because the full action is diffeomorphism invariant, the total stress tensor is covariantly conserved:

$$
\nabla^a\left(T^{(\mu)}_{ab} + T^{(u)}_{ab} + T^{(\mathrm{src})}_{ab} + \sum_i T^{(i)}_{ab} + T^{(\gamma)}_{ab}\right) = 0.
$$

For the active mass-term branch, the exchange terms are explicit:

$$
\nabla^a T^{(\mu)}_{ab}
= \left(-J + \frac12 Z'(\mu)\sum_i m_i^2 \Phi_i^2\right)\nabla_b \mu,
$$

$$
\nabla^a \sum_i T^{(i)}_{ab}
= -\frac12 Z'(\mu)\sum_i m_i^2 \Phi_i^2 \nabla_b \mu,
$$

so the combined metering-plus-massive sector obeys

$$
\nabla^a\left(T^{(\mu)}_{ab} + \sum_i T^{(i)}_{ab}\right) = -J \nabla_b \mu.
$$

A closed source completion therefore requires

$$
\nabla^a T^{(\mathrm{src})}_{ab} = +J \nabla_b \mu
$$

on shell, up to the unit-frame constraint sector. In the explicit source-sector writing above, that condition is no longer an external wish; it is the on-shell exchange law inherited from the same diffeomorphism-invariant action that generates the photon metric, the metering scalar, and the massive-sector mass term. The photon observable is therefore not guessed independently of the massive-sector branch; it is derived from the same metric that solves the sourced Einstein system.

**Propagation speed of $\mu$:** $c_\mu = N(\mu)$. In metered regions ($\mu \gg 1$), $c_\mu \to 1$. In unmetered regions ($\mu \to 0$), $c_\mu \to 0$. The metering density obeys the same causal structure as matter. Horizons are self-reinforcing: once $N = 0$, no $\mu$ dynamics can penetrate.

### Steady State

In the static preferred-frame reduction with the free massive potential $V(\mu) = m^2 \mu^2 / 2$:

$$-\nabla^2 \mu + m^2 \mu = J(x)$$

This is the screened Poisson equation. The screening length $\ell_s = 1/m$ determines how far metering influence extends beyond the source. In the same static frame, the active massive-sector dispersion law reduces to

$$
\omega^2 = k^2 + m_i^2 Z(\mu),
$$

so finite-energy turning occurs at

$$
N(\mu) = \frac{m_i}{E}.
$$

Above threshold, the current workbench uses the corresponding group velocity

$$
v_i(N,E) = \sqrt{1 - \left(\frac{m_i}{E\,N}\right)^2},
$$

the WKB barrier density

$$
\mathcal{I}_i(x;E) = \sqrt{\frac{m_i^2}{N(x)^2} - E^2},
$$

and for a general static path the classical crossing time and excess delay

$$
t_i = \frac{1}{c}\int \frac{dx}{v_i(N(x),E)}, \qquad
\Delta t_i = \frac{1}{c}\int \left(\frac{1}{v_i(N(x),E)} - 1\right)dx,
$$

while for a uniform-lapse segment of physical length $L$ the crossing-time excess is

$$
\Delta t_i(L,N,E) = \frac{L}{c}\left(\frac{1}{v_i(N,E)} - 1\right).
$$

That is the same threshold, barrier, and delay structure used in the current massive-sector and epsilon-regularization calculations.

### Computational Results

Solved on grids of 200 and 1000 points with full nonlinear PDE evolution:

- **Stability:** All linearized eigenvalues $\leq 0$. Zero unstable modes. Horizons are stable.
- **Convergence:** In the tested nonlinear PDE runs, evolution from 5% of steady state converges to within 2.4% of the analytic steady state by $t = 200$.
- **Interior perturbations:** Oscillatory decay (dynamical modes in the metered region).
- **Boundary perturbations:** Frozen. Perturbations at the horizon do not propagate. Eigenfrequencies scale with $N \to 0$.
- **Self-reinforcing horizons:** Signal crossing time through an unmetered gap scales as $\sim \epsilon^{-0.94}$ (theory: $\epsilon^{-1}$). As the gap minimum $N \to 0$, crossing time $\to \infty$.
- **Self-consistent coupling:** Matter-metering feedback converges in ~30--55 iterations in the tested solver. The self-consistent spectrum has higher energies than the prescribed case (narrower effective well), consistent with backreaction.

### Potential Selection

The **free massive potential** $V(\mu) = m^2 \mu^2 / 2$ is the choice adopted in the current model. It produces $\mu \to 0$ in voids (horizons form) and is consistent with the decoherence criterion: if there are no decoherence events ($J = 0$), there should be no metering ($\mu = 0$).

The **symmetry-breaking potential** $V(\mu) = -a\mu^2 + b\mu^4$ has a nonzero vacuum expectation value, giving $\mu > 0$ everywhere including true vacuum. This contradicts the framework's premise, kills horizons, and reduces the bound-state count from hundreds to 2--5.

---

## Cosmological Predictions

### Void Metering and Horizon Formation

A realistic 1D cluster-void-cluster geometry (400 Mpc total, clusters at $\pm 60$ Mpc, void extent ~120 Mpc) was modeled using NFW-like density profiles, temperature-dependent decoherence rates, and the screened Poisson equation.

Key numbers:
- Density contrast: $\rho_{\text{cl}} / \rho_{\text{void}} = 1429$
- Source contrast: $J_{\text{cl}} / J_{\text{void}} = 1.1 \times 10^7$

In the numerical void study, the interior shows a **sharp transition** as a function of the screening mass $m$ and coupling $\alpha$. The controlling parameter is $m \cdot R_{\text{void}}$:

| $m \cdot R_{\text{void}}$ | Behavior |
|---------------------------|----------|
| $< 3$ | $N_{\text{void}} \approx 1$. No observable effect. Metering from clusters fills the void. |
| $3$--$6$ | Transition region. $N_{\text{void}}$ drops from ~1 to ~0 depending on $\alpha$. |
| $> 6$ | $N_{\text{void}} \to 0$. Horizons form. The void interior is causally disconnected for massive particles. |

**Phase diagram** (120-point parameter sweep, $m \in [0.01, 10]$ Mpc$^{-1}$, $\alpha \in [0, 1000]$):

- Horizon-forming ($N_{\text{void}} < 10^{-3}$): **68 / 120** parameter points (57%)
- At the strictest threshold ($N_{\text{void}} < 10^{-6}$): horizons exist for $m \gtrsim 0.36$ Mpc$^{-1}$ at all tested $\alpha$, and never for $m \lesssim 0.17$ Mpc$^{-1}$

Within the studied sweep, the transition is sharp rather than gradual: the tortoise integral numerically separates into horizon-forming and non-horizon-forming regimes.

### Neutrino-Photon Timing Delay

Massive neutrinos ($m_\nu \sim 0.05$--$0.5$ eV) couple to the metering metric. Photons decouple and propagate at $c$. In the same mainline dispersion law, the excess delay for a massive particle traversing a void is

$$
\Delta t_i = \frac{1}{c}\int_{\text{void}} \left(\frac{1}{v_i(N(x),E)} - 1\right) dx,
\qquad
v_i(N,E) = \sqrt{1 - \left(\frac{m_i}{E\,N}\right)^2}.
$$

In the covariant mainline, the sharper statement is that species dependence enters through the turning condition

$$
N(\mu) = \frac{m_i}{E},
$$

so the hierarchy is controlled by the ratio $m_i/E$ rather than by an ad hoc particle label. Heavier or lower-energy massive modes freeze first, lighter or more energetic modes last.

With the current occupancy bridge and the uniform screened reduction $\mu = J/m_\mu^2$, the corresponding activity threshold is explicit:

$$
R_{*,i}
=
\frac{m_\mu^2 L_\perp^2(1+\theta\tau_p)}{\eta_J \tau_p}
\cdot
\frac{1}{\alpha}
\operatorname{arctanh}\!\left(\frac{m_i/E_i - \epsilon}{1-\epsilon}\right),
$$

whenever $\epsilon < m_i/E_i < 1$. If $m_i/E_i \le \epsilon$, then the regularized lapse floor already exceeds the turning threshold and finite-energy exclusion is impossible in that branch. The current workbench now enforces that distinction exactly rather than treating all species with the same qualitative rule.

The same theory layer now carries explicit domain-scaled species benchmarks under the current bridge closure $\Xi_J \approx 1.41 \times 10^{-2}$.

At the current laboratory bound $\alpha = 6.79 \times 10^{-3}$, with $\epsilon = 10^{-2}$, laboratory activity scale $R_{\mathrm{lab}} \sim 10^{40}\ \mathrm{m^{-3}s^{-1}}$, screening length $\ell_s = 0.2\,\mu\mathrm{m}$, and path length `0.2 um`, the uniform benchmark branch gives

$$
N_{\mathrm{lab}} \approx 0.0100000511.
$$

On that benchmark branch:

- a `10 MeV` supernova neutrino with `m = 0.05 eV` is not excluded and picks up only `8.34 x 10^{-29} s` delay across the `0.2 um` gap
- a `1 GeV` electron is not excluded and picks up `8.73 x 10^{-19} s`
- a `1 GeV` proton is excluded, and its turning threshold would require `R_{*,p} \approx 3.33 x 10^{47}\ \mathrm{m^{-3}s^{-1}}`

For a regularized `120 Mpc` reference void with the same $\alpha$, the same closure ratio, $\epsilon = 10^{-2}$, screening length `1 Mpc`, and baryonic baseline activity `R_{\mathrm{void}} \sim 10^{-21}\ \mathrm{m^{-3}s^{-1}}`, the uniform benchmark branch sits at

$$
N_{\mathrm{void}} \approx 0.01.
$$

On that reference branch:

- a `10 MeV` supernova neutrino with `m = 0.05 eV` is not excluded and delays by about `1.55 x 10^3 s`
- a `1 TeV` electron is not excluded and delays by about `1.61 x 10^7 s`
- a `1 GeV` proton is excluded

Those are not full inhomogeneous cosmological transport solutions; they are uniform-domain benchmarks derived from the same occupancy bridge and same-lapse dispersion law. But they already turn the species hierarchy into concrete domain-scale numbers rather than a qualitative slogan.

The same workbench now carries data-anchored messenger windows in [`laboratory_reports/messenger-activity-windows.json`](./laboratory_reports/messenger-activity-windows.json). Three current anchor cases are:

- an `OPERA/CNGS`-like `17 GeV` muon-neutrino messenger over `730 km` with a `10 ns` timing cap
- an `SN1987A`-like `10 MeV` neutrino messenger over `51.4 kpc` with a conservative `3 h` cap
- an `IceCube-170922A`-like `290 TeV` neutrino messenger over `1.75 Gpc` with a loose `1.6 x 10^7 s` association window

Under the present laboratory bound and regularized floor, all three sit in the same regime: the delay-cap-required lapse is already below $\epsilon = 10^{-2}$, so the current branch is trivially safe for these neutrino messengers. Numerically:

- `OPERA/CNGS`-like case: required lapse `1.03 x 10^{-9}`, present lab branch delay `0`
- `SN1987A`-like case: required lapse `7.83 x 10^{-5}`, present reference-void delay `0.661 s`
- `IceCube-170922A`-like case: required lapse `1.29 x 10^{-11}`, present reference-void delay `0`

That is a useful sharpening. Under the current `\alpha` bound and regularized floor, neutrino messengers are not the hard constraint channel. The nontrivial messenger sector remains heavier massive modes, not ultra-relativistic neutrinos.

Selected predictions for a 120 Mpc void:

| $m$ (Mpc$^{-1}$) | $\alpha$ | $N_{\text{void}}$ | Delay |
|-------------------|----------|--------------------|-------|
| 0.05 | 5.0 | 1.000 | 28,001 s |
| 0.10 | 100.0 | 1.000 | 170 yr |
| 0.10 | 20.0 | 0.848 | 8.6 Myr |
| 0.20 | 10.0 | $8 \times 10^{-4}$ | 64 Gyr |
| $\geq 0.50$ | any | $< 10^{-6}$ | effectively infinite |

Over the studied parameter range, the delay appears **binary rather than continuous**: either $N_{\text{void}} \approx 1$ (no interesting delay) or $N_{\text{void}} \ll 1$ (delay exceeds the age of the universe). The detectable window -- delay $> 10$ ns but finite -- is narrow, sitting near the phase boundary.

In this model, that sharpness is a consequence of the phase-transition-like behavior of the void solution.

The epsilon-regularized rerun is now explicit in the current workbench. For the reference same-lapse branch with rest mass `0.1` and conserved energy `8`, increasing $\epsilon$ from `0` to `10^{-2}` drives the WKB action from `2.25544 x 10^5` down to `1.15547 x 10^1` and the tortoise integral from `2.25548 x 10^6` down to `2.53388 x 10^2`, while the bridged 120 Mpc void reference remains below the finite-energy turning threshold (`N_{\text{void}} \approx 0.0100007 < 0.0125` at $\epsilon = 10^{-2}`). So the regularized branch softens the barrier sharply without automatically erasing the horizon-like regime.

The laboratory benchmark now sharpens that statement numerically. If one inserts the current laboratory cutoff $\alpha \le 6.79 \times 10^{-3}$, the benchmark bridge scale $\kappa_J \approx 1.90 \times 10^{-32}\ \mathrm{s\,m^{-2}}$, screening length $\ell_s = 0.2\,\mu\mathrm{m}$, and $\epsilon = 10^{-2}$ into the same-lapse threshold formula above, then the reference mode `(m = 0.1, E = 8)` requires

$$
R_* \approx 4.89 \times 10^{44}\ \mathrm{m}^{-3}\mathrm{s}^{-1}
$$

to reach finite-energy exclusion. So under the present laboratory bound, horizon-scale turning in the regularized branch remains an extreme-activity regime rather than an easy laboratory threshold.

### Differential Gravitational Lensing

In the reduced 1+1D metering geometry, the Ricci scalar $R(x)$ changes sign at the metering boundary:

- **Inside the metered region:** $R > 0$ (convergent lensing)
- **At the boundary:** $R < 0$ (divergent lensing)

This convergence-divergence sign flip is the structural clue carried forward from the reduced model. It is not yet, by itself, the 3+1D photon observable.

The current 3+1D null-geodesic workbench [`three_plus_one_lensing.py`](./three_plus_one_lensing.py) makes the immediate theoretical constraint explicit. Under the mass-term-only photon-decoupling model used in this repository, photons follow the bare metric: the flat decoupled control gives zero direct deflection and zero direct gravitational redshift to numerical precision, while the Schwarzschild control reproduces the weak-field $4M/b$ law at the sub-percent level. A direct metering-induced photon signal therefore requires either an explicit photon-sector completion or a backreaction calculation that feeds the massive-sector metering structure into the metric photons inhabit.

That workbench now contains two explicit nonzero 3+1D completion branches in addition to the flat control.

The first is an Einstein-backreaction branch. A static spherical source $J(r)$ feeds the screened field equation

$$-\nabla^2 \mu + m_\mu^2 \mu = J(r),$$

the resulting profile $\mu(r)$ gives a canonical scalar stress tensor with

$$\rho_\mu = \frac{1}{2}|\partial_r \mu|^2 + \frac{1}{2}m_\mu^2 \mu^2, \qquad p_r = \frac{1}{2}|\partial_r \mu|^2 - \frac{1}{2}m_\mu^2 \mu^2,$$

and the photon metric is then derived from the static spherical Einstein equations

$$m'(r) = 4\pi G r^2 \rho_\mu, \qquad \phi'(r) = \frac{m + 4\pi G r^3 p_r}{r(r - 2m)}.$$

The second is an explicit photon-sector branch, implemented as a disformal-style ansatz

$$g^{(\gamma)}_{\mu\nu} = g^{\text{bare}}_{\mu\nu} - \zeta\,\tanh^2(\alpha_\gamma \mu)\,u_\mu u_\nu + \eta\,\partial_\mu \mu\,\partial_\nu \mu,$$

which in the static spherical workbench becomes a direct photon lapse/radial modification built from the same solved $\mu(r)$ profile.
This is a comparison branch, not yet the adopted mainline completion of the repository.

Under representative weak-field settings, both branches produce nonzero deflection, delay, and redshift, while the flat-decoupled branch remains numerically zero. That gives the project a real 3+1D comparison basis rather than a single placeholder ansatz.

The workbench now also produces axisymmetric lensing profiles and map summaries for cluster-style source families. Representative spherical outputs are in [`theory_maps/einstein-softened-nfw-map.json`](./theory_maps/einstein-softened-nfw-map.json), [`theory_maps/einstein-beta-map.json`](./theory_maps/einstein-beta-map.json), and [`theory_maps/einstein-gaussian-map.json`](./theory_maps/einstein-gaussian-map.json). The spherical mainline result is not the old reduced-model sign flip. For the current spherical Einstein-backreaction branch, cluster-like `softened_nfw` and `beta_model` sources remain convergence-dominated with no resolved negative-$\kappa$ annulus. The Gaussian test can show tiny far-tail negative values after exterior matching, but those sit at the numerical floor of the profile and are not a stable structural sign-flip prediction.

The next step has now been carried out inside the same workbench: multi-component Einstein-branch residual maps with the same `radial_median_bandpass` logic used in the archival program. Representative composite outputs are in [`theory_maps/einstein-single-center-residual-map.json`](./theory_maps/einstein-single-center-residual-map.json), [`theory_maps/einstein-binary-equal-residual-map.json`](./theory_maps/einstein-binary-equal-residual-map.json), [`theory_maps/einstein-bullet-offset-residual-map.json`](./theory_maps/einstein-bullet-offset-residual-map.json), [`theory_maps/einstein-triple-asymmetric-residual-map.json`](./theory_maps/einstein-triple-asymmetric-residual-map.json), and the compact sweep summary [`theory_maps/einstein-composite-residual-sweep.json`](./theory_maps/einstein-composite-residual-sweep.json). That result is sharper. In the viable Einstein branch, the raw convergence map still remains nonnegative across the tested spherical and composite cases. But once the source geometry is genuinely multi-component, the residualized map develops structured positive/negative morphology. The single-center control shows only weak, nearly isotropic filter ringing (`m1 \sim m2 \sim 0.05`), while the equal-binary sweep produces strong quadrupolar residual structure (`m2 \sim 0.36-0.57`) and the bullet-like sweep produces robust dipolar/quadrupolar structure (`m1 \sim 0.21-0.39`, `m2 \sim 0.19-0.44`) across the tested separations and subcluster ratios.

That composite result is no longer only toy. The workbench now has a real-member geometry mode and has been run on the hardened HFF member tables already produced by the archival pipeline. Representative outputs are in [`theory_maps/einstein-abell370-member-geometry-residual-map.json`](./theory_maps/einstein-abell370-member-geometry-residual-map.json), [`theory_maps/einstein-rxcj2248-member-geometry-residual-map.json`](./theory_maps/einstein-rxcj2248-member-geometry-residual-map.json), the cluster summary [`theory_maps/einstein-hff-member-geometry-summary.json`](./theory_maps/einstein-hff-member-geometry-summary.json), and the `top_n` stability check [`theory_maps/einstein-hff-member-geometry-topn-sweep.json`](./theory_maps/einstein-hff-member-geometry-topn-sweep.json). Using the top 25 mass-weighted cluster members from the hardened first-pass tables, both `Abell 370` and `RXC J2248 / Abell S1063` produce sign-changing residual structure while still keeping raw convergence nonnegative. In this first calibration, `Abell 370` gives residual harmonics `m1 \approx 0.19`, `m4 \approx 0.11`; `RXC J2248 / Abell S1063` gives `m1 \approx 0.23`. The `top_n = 10, 25, 40` sweep stays in the same regime for both clusters.

The next bridge is now in place as well: a direct theory-vs-archive comparison mode that rebuilds the HFF residual field from the original `kappa` FITS map using the same `radial_median_bandpass` settings as the archival pipeline, constructs the Einstein member-geometry map on the same member-centroid frame, and compares the two in pixel space and member-score space. Representative outputs are [`theory_maps/abell370-theory-archive-comparison.json`](./theory_maps/abell370-theory-archive-comparison.json), [`theory_maps/rxcj2248-theory-archive-comparison.json`](./theory_maps/rxcj2248-theory-archive-comparison.json), the high-resolution 22-model summary [`theory_maps/hff-model-ensemble-theory-comparison-summary.json`](./theory_maps/hff-model-ensemble-theory-comparison-summary.json), the coarse calibration sweep [`theory_maps/hff-member-geometry-smooth-envelope-calibration-summary.json`](./theory_maps/hff-member-geometry-smooth-envelope-calibration-summary.json), and the rigid-alignment audit [`theory_maps/hff-model-ensemble-alignment-audit-summary.json`](./theory_maps/hff-model-ensemble-alignment-audit-summary.json). That result is more demanding than the morphology-only check, and it has now moved past the first uncalibrated member-only pass. A broad, centroid-contracted smooth envelope built from the same HFF member geometry improves the viable Einstein-branch comparison on the model-family ensemble. At `513^2` map resolution, the 22-model median resolved-sign agreement rises to `\approx 0.497`, the median negative-residual Jaccard overlap rises to `\approx 0.258`, the median member-score Spearman correlation reaches `\approx 0.165`, and the median member pass-agreement reaches `\approx 0.658`. In the primary `cats` maps, the calibrated branch now gives low harmonic mismatch (`L2 \approx 0.062` for `Abell 370`, `0.046` for `RXC J2248 / Abell S1063`) together with stronger structural overlap on the negative side.

The alignment audit sharpens that picture further. Allowing a small rigid registration search (`\pm 20"` in translation, `\pm 20^\circ` in rotation) improves the 22-model median objective from `\approx 1.04` to `\approx 1.30`, the median member-score Spearman from `\approx 0.165` to `\approx 0.449`, the median member pass-agreement from `\approx 0.658` to `\approx 0.712`, the median resolved-sign agreement from `\approx 0.497` to `\approx 0.516`, and the median negative-residual Jaccard from `\approx 0.258` to `\approx 0.275`. But the field-level Pearson correlation still only moves from `\approx 3 \times 10^{-4}` to `\approx 2.5 \times 10^{-3}`. So part of the mismatch is registration-sensitive, but rigid alignment alone does not close the field-level gap.

That changes the theoretical read again. The reduced 1+1 sign flip has still **not** been lifted into a raw 3+1D negative-$\kappa$ prediction of the viable branch. A sign-changing **residual** morphology can now arise within the Einstein-backreaction mainline not only in toy composites and first-pass HFF member geometries, but in a partially calibrated member-plus-envelope geometry that improves the direct archival comparison across both clusters and the 22-model family ensemble. What has **not** happened is a clean pixelwise lock, even after a rigid registration audit. So the immediate theory task is narrower now. It is no longer to ask whether multi-component geometry matters, and no longer simply to add a smooth background. It is to calibrate the viable Einstein branch beyond this first member-plus-envelope stage: amplitude laws, smooth-halo shape, line-of-sight structure, and other cluster-scale degrees of freedom, while also separating physical miscentering from purely comparative registration effects.

The current local 22-family calibration search is now explicit. The first calibrated branch in [`theory_maps/hff-ensemble-calibration-search-summary.json`](./theory_maps/hff-ensemble-calibration-search-summary.json) and the full-resolution rerun in [`theory_maps/hff-ensemble-calibrated-member-envelope-summary.json`](./theory_maps/hff-ensemble-calibrated-member-envelope-summary.json) is a displaced member-plus-envelope Einstein model with `softened_nfw` source profile, smooth-envelope amplitude fraction `0.35`, size multiplier `5.0`, position shrink `0.75`, and a member-frame shift of `10 arcsec` along the principal axis. A subsequent exponent refinement in [`theory_maps/hff-exponent-refinement-summary.json`](./theory_maps/hff-exponent-refinement-summary.json) improves that branch further: the current best member law is amplitude-scaling exponent `0.75` with size-scaling exponent `0.6`. A further smooth-envelope refinement in [`theory_maps/hff-smooth-envelope-refinement-summary.json`](./theory_maps/hff-smooth-envelope-refinement-summary.json) then selects smooth-envelope amplitude fraction `0.25` and size multiplier `6.0` on top of that member law, while a follow-up line-of-sight refinement in [`theory_maps/hff-los-refinement-summary.json`](./theory_maps/hff-los-refinement-summary.json) still leaves `los_amplitude_fraction = 0` as the local winner. A direct smooth-halo offset test in [`theory_maps/hff-smooth-halo-parallel-shift-summary.json`](./theory_maps/hff-smooth-halo-parallel-shift-summary.json) and the higher-resolution check in [`theory_maps/hff-smooth-halo-parallel-shift-highres-comparison.json`](./theory_maps/hff-smooth-halo-parallel-shift-highres-comparison.json) then show that a simple global parallel shift of the smooth envelope is not the missing degree of freedom: the apparent low-resolution `-10 arcsec` improvement does not survive the higher-resolution rerun, and the zero-shift branch remains better.

The next structural probe is now resolved one level deeper. The smooth envelope and LOS layers can now carry an explicitly projected halo shape rather than only isotropic scale factors. In the targeted coarse ensemble check [`theory_maps/hff-elliptic-coarse-summary.json`](./theory_maps/hff-elliptic-coarse-summary.json), the best of the tested branches is an elliptic smooth envelope plus an elliptic shifted LOS component, beating both the isotropic member-envelope baseline and the envelope-only elliptic variants on the 22-model median objective.

The higher-resolution confirmation is now in [`theory_maps/hff-elliptic-best-highres.json`](./theory_maps/hff-elliptic-best-highres.json). Relative to the isotropic displaced member-envelope branch, the current best tested projected branch uses a smooth-envelope axis ratio `q = 0.70`, zero smooth-envelope angle offset from the member principal axis, LOS amplitude fraction `0.10`, LOS axis ratio `q = 0.70`, LOS angle offset `15 deg`, and LOS shift `40 arcsec` along the principal axis. On the aligned 22-family comparison, that branch raises the median objective from `1.4251` to `1.4842`. The strongest gains are structural: aligned resolved-sign agreement improves from `0.4911` to `0.5137`, aligned negative-Jaccard overlap from `0.2811` to `0.3529`, aligned harmonic mismatch drops from `0.1769` to `0.1320`, and aligned member-score Spearman improves from `0.4014` to `0.4141`. The field-level lock is still not there: aligned pixel Pearson stays near zero and does not improve monotonically (`0.0103` to `0.00193`), while aligned pixel Spearman softens slightly (`0.1370` to `0.1266`). So the current best branch is no longer merely anisotropic; it is an elliptic smooth-halo plus shifted elliptic LOS Einstein model, and it improves the overlap structure without yet producing a true pixelwise match.

The current layer-profile split has now been tested explicitly in [`theory_maps/hff-layer-profile-highres.json`](./theory_maps/hff-layer-profile-highres.json). The best tested variant keeps the member layer and smooth envelope on the inherited `softened_nfw` family but switches the shifted elliptic LOS layer to `beta_model`. Relative to the inherited-family baseline, that raises the aligned 22-family median objective from `1.4875` to `1.5221`, the aligned resolved-sign agreement from `0.5157` to `0.5244`, the aligned negative-Jaccard overlap from `0.2973` to `0.3053`, and the aligned member pass-agreement from `0.8750` to `0.8932`. But it is still not the missing final degree of freedom: aligned pixel Pearson remains effectively zero (`2.99 \times 10^{-4}` to `-4.89 \times 10^{-4}`), aligned harmonic mismatch worsens slightly (`0.1303` to `0.1422`), and aligned member-score Spearman softens slightly (`0.4985` to `0.4932`). So the LOS source-family split is a real structural refinement, not a field-level solution.

A direct follow-up sweep over the LOS `beta_model` compactness in [`theory_maps/hff-los-beta-highres.json`](./theory_maps/hff-los-beta-highres.json) does not change that verdict. At scout resolution the LOS branch favored `\beta = 1.0`, but the full-resolution rerun keeps the current `\beta = 0.8` branch as the aligned-objective winner (`1.5221` versus `1.5013`). The higher-$\beta$ variant improves resolved-sign agreement and negative-Jaccard overlap slightly, but not enough to beat the current branch on the full aligned objective. So the remaining gap is not simply a missing LOS `\beta` choice either.

The distributed LOS test is now pinned down as well. The scout pass initially suggested that splitting the LOS layer by `40 arcsec` along the principal axis might help, but the full-resolution confirmation in [`theory_maps/hff-los-distribution-highres.json`](./theory_maps/hff-los-distribution-highres.json) reverses that verdict. The unsplit shifted LOS branch remains the aligned-objective winner (`1.4688` versus `1.4576`), while the split branch improves pixel Spearman and harmonic mismatch slightly but loses on the full aligned objective, negative-side overlap, and member-pass agreement. So distributed LOS structure along that simple principal-axis split is another real probe that does not close the field-level gap.

Simple rigid physical-miscentering by member-frame rotation does not win on its own. In the targeted `513^2` confirmation [`theory_maps/hff-member-rotation-highres.json`](./theory_maps/hff-member-rotation-highres.json), the zero-rotation branch remains the aligned-objective winner (`1.5221`) over both `+5 deg` (`1.5001`) and `-5 deg` (`1.4079`). The `+5 deg` branch does improve member-score Pearson (`0.2960` to `0.4013`) and member-score Spearman (`0.4932` to `0.5053`), but it still loses on the full aligned objective. So rigid rotation alone is not the missing degree of freedom.

The same is true of the richer smooth-halo family in isolation. At `257^2`, replacing the inherited smooth envelope with a `beta_model` smooth halo sharpened the structural overlap. But the `513^2` confirmation in [`theory_maps/hff-smooth-source-family-highres.json`](./theory_maps/hff-smooth-source-family-highres.json) returns the inherited smooth-envelope family as the aligned-objective winner (`1.5221` versus `1.4997`). The `beta_model` smooth halo improves aligned resolved-sign agreement (`0.5244` to `0.5694`), negative-Jaccard overlap (`0.3053` to `0.3486`), and harmonic mismatch (`0.1422` down to `0.1141`), but it gives up member-score Spearman (`0.4932` to `0.4236`) and member-pass agreement (`0.8932` to `0.8717`).

The first branch that actually survives those pieces in combination is now explicit in [`theory_maps/hff-targeted-variant-highres.json`](./theory_maps/hff-targeted-variant-highres.json). In that targeted `513^2` comparison, the aligned-objective winner is a `+5 deg` member-frame rotation together with a `beta_model` smooth envelope (`\beta = 0.8`), keeping the shifted elliptic LOS `beta_model` branch in place. Relative to the old baseline branch, the aligned 22-family medians move as follows:

- objective: `1.5221` to `1.5706`
- pixel Pearson: `-4.89 x 10^-4` to `3.12 x 10^-4`
- pixel Spearman: `0.1199` to `0.1399`
- resolved-sign agreement: `0.5244` to `0.5601`
- negative-Jaccard overlap: `0.3053` to `0.3427`
- harmonic `L2`: `0.1422` down to `0.1149`
- member-score Pearson: `0.2960` to `0.3593`
- member-score Spearman: `0.4932` to `0.5142`
- member-pass agreement: `0.8932` to `0.8894`

So the lensing mainline moved again. The missing degree of freedom is not simple LOS splitting, not simple rotation alone, and not smooth-halo family alone. But the combination of a mild physical miscentering rotation with the richer smooth-halo family does produce a stronger field-level comparison. Pixelwise lock is still weak, yet the current best branch now improves both overlap structure and member-level ordering without sacrificing the full aligned objective.

The next structural probe is now pinned down too. A broader harmonic-modulated halo search in [`theory_maps/hff-harmonic-variant-highres.json`](./theory_maps/hff-harmonic-variant-highres.json) found a non-mainline branch with smooth-envelope quadrupole modulation `m2 = 0.3` and phase offset `15 deg` that can raise aligned pixel Pearson to about `2.27 x 10^{-3}`, but it loses sharply on the full aligned objective (`1.4245`) and on harmonic mismatch. A direct follow-up on top of the current best branch in [`theory_maps/hff-targeted-harmonic-midres.json`](./theory_maps/hff-targeted-harmonic-midres.json) keeps the baseline winner: adding smooth-envelope quadrupole modulation with `m2 = 0.15` or `0.30` lowers the aligned median objective from `1.4927` to `1.4283` or `1.4342`. So simple low-order angular modulation is not the missing field-lock degree of freedom either.

The next exact comparison has now closed the transverse-miscentering question one level further. In [`theory_maps/hff-targeted-miscentering-highres.json`](./theory_maps/hff-targeted-miscentering-highres.json), the live branch was tested directly against nearby challengers built from physical member-frame and large-scale transverse offsets: `\pm 10` arcsec member-frame perpendicular shifts, a `20` arcsec member-frame parallel shift, `10` arcsec smooth-envelope perpendicular and parallel shifts, a `10` arcsec LOS perpendicular shift, a `20` arcsec perpendicular LOS split, and an anisotropic smooth-envelope shrink variant. At mid resolution some of those variants looked competitive enough to promote, but the `513^2` confirmation keeps the current baseline winner. The strongest challengers slightly improve narrow diagnostics such as aligned pixel Pearson, aligned pixel Spearman, or aligned member-score Spearman, but they lose on the full aligned objective. Numerically, the baseline branch remains ahead with aligned median objective `1.5138`, against `1.5123` for a `+10` arcsec smooth-envelope perpendicular shift, `1.5086` for a `+10` arcsec LOS perpendicular shift, and `1.4802` for a `20` arcsec member-frame parallel shift. So simple transverse physical miscentering is not the missing field-lock degree of freedom either.

A direct global-strength sweep in [`theory_maps/hff-backreaction-scale-summary.json`](./theory_maps/hff-backreaction-scale-summary.json) then shows that the current reference normalization is already near the local optimum of the present scaffold: weakening or strengthening the overall Einstein backreaction does not open a cleaner field-level fit. More sharply still, the factorization audit in [`theory_maps/hff-backreaction-factorization-invariance.json`](./theory_maps/hff-backreaction-factorization-invariance.json) shows exact invariance of the 22-family comparison when the composite scale

$$
S_{\mathrm{eff}} = \texttt{gravity\_scale}\,\texttt{density\_scale}\,\texttt{source\_amplitude}^2
$$

is held fixed. So within the current Einstein scaffold, the lensing calibration determines a composite backreaction scale, not the separate source amplitude and stress-energy normalizations.

The calibration tool now accepts bridge-driven source amplitudes as well: the Einstein search can take either an absolute bridge conversion $\kappa_J$ plus target activity density $R$, or a reference source amplitude together with a reference-to-target activity ratio. That makes the lensing scaffold compatible with the same occupancy law used in the laboratory and void sections, even though the present comparison still constrains only the composite backreaction scale rather than the bridge parameters separately.

So the viable branch is now materially better calibrated than the first member-plus-envelope pass. The remaining gap is not missing morphology, not missing halo-shape freedom, not missing LOS freedom, not missing a trivial global rescale, and not simply missing a LOS source-family split. The remaining gap is a real field-level one: the current Einstein branch is now competitive on sign structure, overlap structure, and member-level ordering, but not yet on full pixelwise lock.

**Observational target:** Galaxy cluster lensing residuals. The archival program tests whether the sign-flip-like morphology survives model-family variation and whether any decoherence-relevant proxy carries conditional signal beyond stellar mass and local convergence. That is a morphology-and-confounder program, not yet a completed derivation of a direct photon-sector effect.

The current archival program works on the two public HFF cluster fields (`Abell 370`, `RXC J2248 / Abell S1063`) using the official HFF property catalogs, BUFFALO v2 photometric proxy catalogs, public MUSE Lensing Clusters spectroscopic catalogs, a hardened `radial_median_bandpass` residual construction, deeper external spectroscopy for membership hardening, public `MUSE-DEEP` core-cube extractions, and 22 independent frontier-model `kappa` maps. The sign-flip-like morphology is reproducible across model families, so the observational question is no longer whether the morphology exists in the archive; it is which explanatory axis survives the controls and replicates across clusters.

At present, `local_kappa` is the strongest and most stable predictor. Across the 22 model-family runs it is positive in `21/22`, significant in `19/22`, and has median partial correlation `0.3604`. Stellar mass remains a real axis as well: positive in `19/22`, significant in `16/22`, with median partial correlation `0.3740`. In the primary two-cluster pass, both combine at `p \approx 5.35 \times 10^{-5}`.

The BUFFALO-based continuum proxy does not survive the hardened residual. Across the 22 model-family runs it is positive in `8/22`, significant in `1/22`, with median partial correlation `-0.0316`.

The stronger spectroscopic pass is more informative, but it does not yet produce replication. MUSE-derived line-count and optical-complexity proxies are positive in `Abell 370` and negative in `RXC J2248 / Abell S1063` in the primary two-cluster pass. Across the 22 model-family runs, line count is positive in `14/22` and significant in `7/22`, with median partial correlation `0.1002`; optical complexity is positive in `14/22` and significant in `2/22`, with median partial correlation `0.1040`; emission strength is positive in `15/22` and significant in `2/22`, with median partial correlation `0.0786`. Those are real non-null structures, but not same-sign replication across the two clusters.

The first direct strong-line metallicity pass has now also been carried through using public MUSE line catalogs. In `Abell 370`, `N2` and Marino et al. (2013) oxygen-abundance estimators are mildly positive but not significant once mass, redshift, magnification, radius, and `local_kappa` are controlled. In `RXC J2248 / Abell S1063`, the same estimators are effectively sample-limited under the current public coverage. Across the 22 model-family runs, `muse_n2` is positive in `6/22` and significant in `0/22` with median partial correlation `0.0304`; the Marino `O3N2` abundance estimator is positive in `7/22` and significant in `0/22` with median partial correlation `0.0672`; the Marino `N2` abundance estimator is positive in `4/22` and significant in `0/22` with median partial correlation `-0.0480`.

The public matched-depth core-cube route has now been tested as well. Secure-redshift aperture extractions from the public `MUSE-DEEP` core cubes do not rescue the metallicity lane: the cubes cut off near `8850-8900 A`, so `Halpha/[N II]` fall out of band or effectively out of band at the cluster redshifts, and the remaining blue-line support is too thin to carry a two-cluster abundance test. In the current deep-cube extraction cache, `Abell 370` yields only `2` usable `R23` objects and `9` usable `O32` objects, while `RXC J2248 / Abell S1063` yields `0` usable `R23` objects and `1` usable `O32` object; `N2`, `O3N2`, and Marino-style abundance proxies remain unsupported in both clusters. That closes the observational proxy lane on current public data. Two later observational pursuits remain: redder symmetric spectroscopy across both clusters with usable `Halpha`, `[N II]`, `Hbeta`, and `[O III]`, and direct per-galaxy abundance tables mapped one-to-one into the HFF member catalogs. The active center of effort now shifts to the 3+1D lensing derivation.

### GW170817 Compatibility

GW170817/GRB 170817A constrains $|c_{\text{GW}} - c_{\text{photon}}|/c < 10^{-15}$.

Since the metering coupling enters through the mass term, all massless fields (photons, gravitons, gluons) decouple. Gravitational waves propagate at $c$ everywhere, regardless of $\mu$. Within this mass-term-only coupling model, the GW170817 constraint is satisfied for any $(m, \alpha)$, with no parameter tuning required. The explicit photon-sector comparison branch is retained only as a constrained counterexample and not as the current mainline.

---

## Observational Proposals

### A. Zeno Spectral Test (Completed)

A quantum system under spatially varying measurement strength exhibits discrete bound-state spectral structure matching the $1/N^2$ potential, rather than generic Zeno freezing. Realizable with current ion trap or superconducting qubit technology. **Three IBM Quantum hardware runs completed** (Section: [IBM Quantum Hardware](#ibm-quantum-hardware)).

### B. Anomalous Casimir Effect

The WEC violation predicts an additional negative-energy contribution in regions of high metering density. That lane is no longer just a proposal. A first executable laboratory workbench now exists in [`laboratory_bounds.py`](./laboratory_bounds.py), with canonical reports in [`laboratory_reports/casimir-benchmark-report.json`](./laboratory_reports/casimir-benchmark-report.json), [`laboratory_reports/casimir-differential-design.json`](./laboratory_reports/casimir-differential-design.json), [`laboratory_reports/casimir-planar-slab-geometry.json`](./laboratory_reports/casimir-planar-slab-geometry.json), and [`laboratory_reports/casimir-source-bridge-report.json`](./laboratory_reports/casimir-source-bridge-report.json).

Under the repository's current benchmark model:

- the strongest executed Casimir pressure null gives $\alpha \leq 6.79 \times 10^{-3}$
- a signal at $\alpha = 10^{-2}$ would require an active-vs-inert pressure difference of about `1.30 mPa` at `0.2 um`
- a signal at $\alpha = 10^{-3}$ would require about `13.0 uPa` at the same separation

The first reference differential design is now explicit as well: `13.0 pN` at `0.2 um` for a `1 mm^2` active area when $\alpha = 10^{-3}`, together with common-mode rejection `6.25 x 10^4` and uncompensated gap-drift tolerance `0.8 pm`. Relative to the `0.3 pN` static-force absolute error quoted in the Decca 2003 Casimir measurement, that reference target is larger by a factor `\approx 43.3`.

The workbench now also carries a first geometry-dependent planar slab model in [`laboratory_reports/casimir-planar-slab-geometry.json`](./laboratory_reports/casimir-planar-slab-geometry.json). In that model, two active `0.05 um` source slabs feed the screened field equation with screening length `0.2 um`, and the slab interaction law is derived in closed form as

$$P_{\text{slab}}(g) = \frac{J_0^2}{2m_\mu^2} e^{-m_\mu g} \left(1 - e^{-m_\mu t}\right)^2.$$

The finite-domain numerical solver in the workbench reproduces that closed form at the `10^{-5}` level across the current gap scan, so the planar interaction shape is now analytically pinned down inside the linear screened model. The main consequence is that the gap dependence is no longer forced to track pure Casimir `1/a^4` scaling. Under this calibrated slab model:

- at `0.3 um`, the geometry modifier is `3.07` relative to the pure benchmark extrapolation
- at `0.5 um`, the geometry modifier is `8.72`

The remaining laboratory unknown is no longer hidden in a generic calibration constant. For the current default slab geometry, matching the `0.2 um` reference amplitude requires an effective slab source strength `J_0 \approx 1.90 x 10^8` in the screened-field equation, which corresponds to a bridge conversion

$$\kappa_J \approx 1.90 \times 10^{-32}\ \text{s m}^{-2}$$

relative to the repository's laboratory decoherence-activity scale `R_lab ~ 10^40`. Relative to the naive causal screening scale `1/(c \ell_s)`, that required conversion is smaller by a factor `\approx 1.14 x 10^{-30}`. So the naive source interpretation is not merely incomplete; in this form it is effectively ruled out.

The new source-bridge report shows that the exact slab law already selects a much sharper alternative. If the static source is written as

$$J(x) = \eta_J \frac{\tau_p}{L_\perp^2} R(x),$$

then a Planck-tick occupancy bridge with $\tau_p = t_P$ gives two clean reference statements:

- if $L_\perp = \ell_s = 0.2\,\mu\text{m}$, the required bridge efficiency is $\eta_J \approx 1.41 \times 10^{-2}$
- if $\eta_J = 1$, the required transverse coarse-graining scale is $L_\perp \approx 1.68\,\mu\text{m}$

That does not finish the covariant source theory, but it does move the bridge from an impossible suppression problem to a narrow and structurally specific occupancy family.

The lensing side now sharpens what that bridge can and cannot fix. Within the current Einstein-backreaction scaffold, the 22-family HFF comparison is exactly invariant across the tested refactorings

$$
(\texttt{source\_amplitude}, \texttt{density\_scale}) =
(5\times 10^{-4}, 4\times 10^{-6}),
(10^{-3}, 10^{-6}),
(2\times 10^{-3}, 2.5\times 10^{-7}),
$$

all at fixed

$$
S_{\mathrm{eff}} = \texttt{gravity\_scale}\,\texttt{density\_scale}\,\texttt{source\_amplitude}^2 = 10^{-12}.
$$

So the present lensing branch does not separately determine the source amplitude and stress-energy normalization. It fixes only the composite backreaction scale until the covariant source law is closed.

So this remains the shortest current path to a hard physical bound, but it is no longer vague. The next laboratory move is a purpose-built active-vs-inert Casimir differential with `O(10 uPa)` pressure reach and `O(10 pN)` force reach near `0.2 um`, together with wider-gap follow-up runs where the exact screened-slab law predicts a slower falloff than the raw `1/a^4` benchmark guess, and a covariant fixation of the occupancy bridge parameters $(\eta_J, \tau_p, L_\perp)$ that turns the current bridge family into a true physical prediction.

### C. Cosmological Metering Signature

The effective lapse function at early times (low metering density) differs from today. This produces a scale-dependent modification to the CMB power spectrum correlated with the information-processing capacity of matter at each epoch.

### D. Geometric Observer Detection

If a consistent 3+1D completion yields a nontrivial photon observable, regions containing observers could imprint measurably different physical properties:

- **Anomalous gravitational lensing:** a metered region could lens light differently than an equal-mass control once the photon sector or backreaction is specified.
- **Redshift anomaly:** a metering-dependent frequency shift would have to be derived in the metric photons actually inhabit.
- **Thermal halo:** The analog Unruh temperature at the metering boundary peaks at a characteristic value set by the metering gradient.
- **Correlation signature:** A true metering signal shows correlated anomalies across all three channels (lensing, redshift, thermal). False positives from unrelated astrophysics do not produce this correlation.

If any of these channels proved real, they would suggest a new SETI-style methodology: searching for observers by their geometric imprint on spacetime, not by engineering artifacts or communication attempts.

### E. Void Lensing Residuals (Archival)

The convergence-divergence sign-flip morphology at void boundaries is proposed here as a structural signature that differs from standard CDM subhalo convergence patterns and from simple amplitude-only modified-gravity distortions. Test: lensing residuals in existing cluster surveys, correlated with proxies for information-processing density. Requires no new observations.

The working archival lensing program now has three pieces:

1. establish the morphology itself under model-family variation,
2. force every proposed proxy to compete against mass and local convergence under explicit controls,
3. only expand the cluster sample after a proxy survives that conditional test.

Within that program, the hardened residual baseline is now in place, stronger public spectroscopic catalogs have been tested, a first strong-line metallicity pass has been completed where line coverage allows, and the public matched-depth core-cube route has been pushed through secure-redshift aperture extraction. On current public data, that observational proxy lane is exhausted. The later observational pursuits are redder symmetric spectroscopy across both clusters and direct abundance products with one-to-one HFF mapping. The active lensing task is now the 3+1D derivation of the observable itself.

---

## Open Questions

### 3+1D Lensing Observable

The project now needs the actual lensing observable derived from a consistent 3+1D completion, not further proxy reshuffling on the current public archive. The central task is to write the covariant theory cleanly enough that photon propagation, lensing, and any claimed redshift effect can be computed from the field equations rather than inferred from lower-dimensional geometric intuition.

A concrete workbench for that step now exists in [`three_plus_one_lensing.py`](./three_plus_one_lensing.py). It evaluates null geodesics for static spherical photon metrics, checks the flat-decoupled control, benchmarks against Schwarzschild, solves a sourced spherical $\mu(r)$ profile from the screened field equation, derives an Einstein-backreaction branch from the resulting stress tensor, and evaluates an explicit photon-sector coupling branch on the same profile. At present that workbench says three sharp things:

1. flat photon decoupling yields no direct metering-induced photon deflection or redshift;
2. a field-equation-derived Einstein-backreaction branch yields a computable weak-field photon signal;
3. an explicit photon-sector coupling branch yields a separate computable nonzero photon signal on the same source profile.

The workbench now also carries a first constraint audit against GW170817-style speed bounds, solar-system light-bending and Shapiro-delay $\gamma$ bounds, and multimessenger weak-equivalence-principle bounds. The result is asymmetric. The explicit photon-sector branch is easy to generate but hard to keep: under the representative weak-field profile now in the repository, the default direct-coupling choice fails the VLBA, Cassini, and GW170817 $\Delta\gamma$ tests by orders of magnitude. Surviving parameter space is not broad; it lies on a narrow tuning line relating the temporal and radial photon couplings. In the current audit at reference radius $r = 360$ for the default sourced profile, a temporal coupling of `0.02` requires a radial coupling near `35.95` to land at $\gamma \approx 1$, while the current illustrative value `0.1` is far outside that allowed band.

GW170817 is also sharp about asymptotic background coupling. For a homogeneous nonzero metering background $\mu_\infty$, the direct photon branch picks up an asymptotic speed shift unless the temporal coupling is extremely small. In the current audit code this appears as

$$\zeta_{\max} \sim \frac{2\,\delta c_{\max}}{\tanh^2(\alpha_\gamma \mu_\infty)},$$

with $\delta c_{\max} \sim 3 \times 10^{-15}$. Even a tiny background such as $\mu_\infty = 10^{-3}$ drives the allowed temporal coupling down to order `6e-9`.

The Einstein-backreaction branch is in better shape. In the exterior it is a GR-derived branch and relaxes toward Schwarzschild with $\gamma = 1$. That makes it the current mainline candidate. In spherical cluster-style maps, the branch remains convergence-dominated and does not produce a resolved raw sign-flip morphology. In multi-component residual maps it does produce the right sign-changing residual class while keeping the raw convergence map nonnegative. On real HFF member geometry, a displaced member-plus-envelope model improved the 22-family comparison; the next refinement showed that projected halo shape matters; and the current best tested structural branch is an elliptic smooth envelope plus a shifted elliptic LOS component. In the current high-resolution comparison that branch raises the aligned 22-family median objective to `1.4842`, improves resolved-sign agreement to `0.5137`, raises negative-Jaccard overlap to `0.3529`, and lowers harmonic mismatch to `0.1320`. But the field-level gap remains: pixelwise correlation still stays near zero, and the factorization audit still shows that the present lensing comparison fixes only the composite scale $\texttt{gravity\_scale}\,\texttt{density\_scale}\,\texttt{source\_amplitude}^2$, not the source and stress-energy normalizations separately. So the current open problem is no longer whether the viable branch can make the residual class. It can. The open problem is driving that branch to a quantitatively competitive field match with a fully closed source law.

### The Coupling Constant

$\alpha$ is the central unknown. The hardware fit ($\alpha \approx 1.16$) applies to the lattice model, not to spacetime. The internal WEC benchmark ($\alpha < 0.133$) is a consistency cutoff, not a laboratory result. The repository now also contains a first executed measurement-anchored Casimir null bound, with strongest current benchmark value $\alpha \leq 6.79 \times 10^{-3}$ under the present anomaly model, plus a first explicit active-vs-inert differential design target at the `10^{-3}` level, a first exact planar screened-slab interaction law for the laboratory geometry, and a first viable occupancy-style source bridge family. The dimensional estimate ($\alpha_{\text{SI}} \sim \ell_P^3$) still suggests the physical coupling may be extraordinarily small. All cosmological predictions depend on $\alpha$, and a direct first-principles physical extraction of $\alpha$ now reduces, on the laboratory side, to fixing the covariant law for $(\eta_J, \tau_p, L_\perp)$ rather than to inventing a bridge from scratch.

### Asymptotic Geometry

The Ricci scalar grows as $R \sim -x^2 / (4\sigma^2)$ in the unmetered region and diverges. The horizon is profile-independent, but the asymptotic curvature is unbounded. This may signal that the unmetered region is not a valid spacetime at all -- not just "no duration" but "no geometry." The curvature diverges because the description is being pushed past its domain of validity.

### Relationship to the Measurement Problem

The framework distinguishes "measured" from "unmeasured" at a geometric level. This is closer to objective collapse theories (GRW, Penrose) than to many-worlds or Copenhagen. In many-worlds, every branch has observers relative to itself; the framework would require reformulation. The framework is not interpretation-neutral.

### Information-Theoretic Content

The bound-state count (~69 radial in 1+1D, ~600 at 1000 sites, ~17 bits in 3+1D with angular modes) suggests connections to holography:

- Does the bit count scale with boundary area (holographic) or volume?
- Does it match Bekenstein or Bousso bounds?
- What information is geometrically confined?

### Energy Dependence

The mass-term coupling means $m_{\text{eff}} = m_0 / N(x)$. Heavier particles have larger absolute mass divergence in voids. The species hierarchy for void propagation is: heavy particles freeze first, light particles last. Neutrinos (lightest massive) are least affected. The current same-lapse workbench now enforces the regularized-floor version of that statement exactly: if a species threshold $m_i/E_i$ drops below the lapse floor $\epsilon$, exclusion is impossible in that regularized branch; otherwise the required activity threshold follows from the explicit $R_{*,i}$ formula above.

### Deferred Covariant Branch: $N^{-4}$ Massive-Sector Completion

A distinct theoretical branch has now been identified but is being kept separate from the active line of work. In that branch, the massive sector couples to a static metering metric

$$ds_m^2 = -N(x)^2 dt^2 + h_{ij}(x) dx^i dx^j$$

with a local mass-shell term scaled as

$$\mathcal{L}_{\text{massive}} \sim -\frac{1}{2}\sqrt{-g_m}\left(g_m^{\mu\nu}\partial_\mu \Phi \, \partial_\nu \Phi + \frac{m_0^2}{N(x)^4}\Phi^2\right)$$

so that

$$m_{\text{eff}}(x) = \frac{m_0}{N(x)^2}.$$

This branch matters because if the same $N$ serves both as lapse and as the field controlling the mass term, then the simpler $m_0^2/N^2$ scaling does not by itself create a horizon-driven forbidden region in the local Hamilton-Jacobi equation: the mass term and the redshifted energy term scale the same way. By contrast, the $N^{-4}$ branch gives a finite-energy turning surface

$$N_*(E) = \frac{m_0}{E},$$

so for $N < N_*(E)$ the massive mode is classically excluded. Near a simple zero of the lapse, the WKB barrier grows rapidly and the near-horizon wavefunction is exponentially suppressed. In that sense, the $N^{-4}$ branch gives a cleaner local covariant realization of horizon-adjacent confinement than the current same-lapse $1/N^2$ story.

This is **not** the active branch of the repository. The current theorem spine, lattice numerics, hardware comparison, and README claims are organized around the tested $1/N^2$ potential. Moving to the $N^{-4}$ branch would therefore require a deliberate reworking of the covariant completion, the continuum spectral analysis, and the interpretation of the existing lattice evidence. It is being recorded here as a deferred branch to evaluate only after the present mainline program has been carried through.

---

## Negative Results

These are documented for completeness and to prevent others from repeating dead ends.

**Geodesic normalization (v1):** The lapse-function approach to conditional proper time fails. Geodesic equations enforce $d\tau/d\lambda = 1$ automatically. Worldlines accumulate proper time in unmetered regions regardless of $N$. The mass-gap formulation (divergent effective mass) is the correct mechanism.

**Entropy test bug (v2):** The entropy-based duration test returned zero transitions because the random walk preserved the distribution used in the entropy calculation. Fixed by using a sliding window over visit history.

**Born-Oppenheimer decomposition:** The BO decomposition of the Wheeler-DeWitt solution never cleanly demonstrated emergent time. The coefficient of variation remained above 0.85 in all tests. The toy minisuperspace model (120 x 120 grid) is too crude. This does not affect the main results (which are independent of the BO analysis).

**Public strong-line and matched-depth metallicity passes:** The current public MUSE line catalogs are sufficient to construct a first `N2` / Marino-style abundance pass in `Abell 370`, but they remain too sparse in `RXC J2248 / Abell S1063` to support a symmetric two-cluster test. The public `MUSE-DEEP` core cubes are also now in hand and scriptable, but at the cluster redshifts their wavelength coverage ends before a viable two-cluster `Halpha/[N II]` abundance test, and the surviving blue-line support is too thin to substitute. The result is informative as a data-adequacy and wavelength-coverage bound. The current public-data proxy lane is therefore exhausted; later observational pursuit requires redder symmetric spectroscopy or direct abundance tables.

---

## Computational Verification

The computational components of the project are checked using independent tools:

- **NumPy/SciPy:** Wheeler-DeWitt eigensolve, geodesic integration, density matrix evolution, PDE integration
- **Wolfram Mathematica 14.3:** 23 symbolic verification modules covering Lagrangian analysis, Christoffel symbols, Ricci scalar, tortoise coordinates, WKB tunneling, bound states, energy conditions, stability analysis, 3+1D extension
- **Qiskit 2.3.0 + IBM Quantum:** 3 hardware runs on `ibm_torino` (2026-02-23)

Simulator results span two orders of magnitude in lattice resolution (6 to 1000 sites) with spectral accuracy improving from 7% to 0.06%, confirming convergence to the continuum limit.

---

## Status

| Component | Status |
|-----------|--------|
| Core theorems (horizon, confinement, proper time) | Abstract theorem spine formalized in Coq; dimensional cases computationally checked |
| Zeno spectral test (simulator, 6--1000 sites) | Complete, 0.06% accuracy |
| Zeno spectral test (IBM hardware, 3 runs) | Complete, revival observed in the tested lattice model |
| Potential shape discrimination | 6-site hardware comparison complete; $1/N^2$ gives the best fit among tested shapes |
| Metering field equation (PDE + stability) | Numerical study complete; stable in the tested setups |
| Photon decoupling | Model consequence of the mass-term-only coupling |
| 3+1D lensing observable | Flat control, Einstein-backreaction branch, and explicit photon-sector branch are computable; the direct photon branch is heavily constrained; the Einstein branch is the current mainline candidate; spherical cluster-style maps remain convergence-dominated; multi-component, HFF-member, and calibrated member-plus-envelope Einstein maps generate sign-changing residual morphology without raw negative convergence; the current best tested structural branch is a `+5 deg` member-frame rotation plus a `beta_model` smooth envelope (`\beta = 0.8`) together with the shifted elliptic LOS `beta_model` branch; simple principal-axis LOS splitting, simple rigid rotation alone, smooth-`beta_model` replacement alone, simple low-order harmonic halo modulation, and simple transverse miscentering probes do not beat that combined branch; pixelwise field matching remains weak |
| Decoherence rate source criterion | Dynamic activity density $R = \sum n_i \gamma_{D,i}$ adopted; the static screened source is now written through an occupancy transport law $\nabla_a(Q u^a) = R - Q/\tau_p$ with $J = \eta_J Q/L_\perp^2$; the theory-native closure reduces this to the dimensionless ratio $\Xi_J = \eta_J \chi_\tau/\chi_L^2$ once $(\tau_p,L_\perp)=(\chi_\tau t_P,\chi_L \ell_s)$ is imposed; in the minimal natural reduction, $\Xi_J$ is the screening-cell occupancy fraction, with current benchmark value `1.41 x 10^-2`; exact slab work fixes the benchmark bridge scale; the current reference Einstein source amplitude `10^{-3}` corresponds to bridge-side activity scale `R \approx 5.26 \times 10^{28}\ \mathrm{m}^{-3}\mathrm{s}^{-1}` at that benchmark bridge; the new cross-domain audit rules out a universal literal one-cell Planck-tick closure, so the surviving bridge problem is a renormalized effective reduction rather than a raw microscopic cell count; the current Einstein lensing scaffold still fixes only the composite backreaction scale $\texttt{gravity\_scale}\,\texttt{density\_scale}\,\texttt{source\_amplitude}^2`, but the calibration tool now accepts bridge-derived source amplitudes from either absolute bridge conversion or reference-to-target activity ratios |
| Void metering + cosmological predictions | Numerical parameter study complete; observational status open |
| Archival lensing program (2 clusters, 22 model families, BUFFALO + MUSE + deeper external spectroscopy + public MUSE-DEEP core cubes) | Morphology replicated under the hardened residual; local convergence and mass remain stable; the public-data proxy lane is exhausted without a replicated abundance-style signal; active focus shifts to the 3+1D lensing derivation |
| Coq formalization | Abstract theorem spine compiles; admissible-path and integral bounds, massive-threshold results, mass- and energy-hierarchy monotonicity, and local/interval regularized-floor theorems are formalized in the monolithic file |
| Physical $\alpha$ laboratory program | Casimir benchmark null executed; strongest current benchmark cutoff $\alpha \leq 6.79 \times 10^{-3}$; active-vs-inert differential target, exact planar slab law, and first viable occupancy bridge family defined; remaining blocker is fixing $(\eta_J, \tau_p, L_\perp)$ covariantly; direct first-principles measurement remains open |

---

## License

MIT
