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
- **Measurement-anchored laboratory benchmark:** [`laboratory_bounds.py`](./laboratory_bounds.py) now maps the repository's current Casimir-scale anomaly benchmark onto published precision Casimir-pressure nulls. Under the present benchmark model, the strongest executed laboratory null in the repository is a Casimir bound of $\alpha \leq 6.79 \times 10^{-3}$, with a canonical report in [`laboratory_reports/casimir-benchmark-report.json`](./laboratory_reports/casimir-benchmark-report.json).
- **Archival observational program:** the cluster-lensing section currently runs on two HFF fields, a hardened radial-median-bandpass residual, BUFFALO v2 photometric catalogs, public MUSE spectroscopic catalogs, deeper external spectroscopy for membership hardening, public `MUSE-DEEP` core-cube extractions, and 22 independent lens-model families. The sign-flip morphology survives those checks. The most stable explanatory axes are local convergence and stellar mass. Generic spectroscopic structure, line-ratio proxies, first strong-line metallicity estimators, and the current public matched-depth core-cube route do not yield a same-sign replicated abundance-style proxy across both clusters under the present controls.
- **Open empirical question:** the cosmological interpretation, observer-detection ideas, and any claim that the measured morphology tracks decoherence-rate proxies beyond mass and local convergence remain open. The observational proxy lane is exhausted on current public data. Two later observational pursuits remain well-defined: redder symmetric spectroscopy for both clusters with usable `Halpha`, `[N II]`, `Hbeta`, and `[O III]` coverage, and direct per-galaxy abundance tables with one-to-one mapping into the HFF member catalogs. Until such data exist, the active focus shifts to deriving the 3+1D lensing observable from the field equations.

---

## Claim Ladder

The repository now supports a cleaner claim hierarchy:

- **Established within the current project:** the mass-gap mechanism is the working core. The abstract theorem spine in [`meteringmetric.v`](./meteringmetric.v) supports horizon formation, confinement-barrier structure, and proper-time vanishing in the present theorem-level scope. The lattice and hardware program supports discrete bound-state structure for the tested $1/N^2$ potential. The archival lensing program supports a reproducible sign-flip-like residual morphology, with local convergence and stellar mass as the stable explanatory axes on current public data.
- **Supported but not finished:** a viable 3+1D Einstein-backreaction branch exists in [`three_plus_one_lensing.py`](./three_plus_one_lensing.py). In that branch, raw convergence remains nonnegative in the tested cluster-style maps, but sign-changing residual morphology appears in multi-component, HFF-member, and calibrated member-plus-envelope geometries. Direct archive comparison improves under that calibration, and a rigid-alignment audit shows that part of the remaining mismatch is registration-sensitive.
- **Open:** a full 3+1D observable derivation is still needed. The current theory program has not yet turned the viable branch into a quantitatively competitive field-level match to the archival residual maps. A continuum finite-dimensionality theorem remains open. A direct first-principles physical measurement of the spacetime coupling $\alpha$ remains open, even though the repository now contains a measurement-anchored Casimir benchmark null bound. On the observational side, the next clean empirical step would require redder symmetric spectroscopy or direct abundance tables for both clusters.
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

### Recommended: Decoherence Rate Density

$$J(x) = \alpha_J \sum_{\text{massive } i} n_i(x) \cdot \gamma_{D,i}(x)$$

where $n_i$ is the number density of massive species $i$ and $\gamma_{D,i}$ is the local decoherence rate.

**Why this criterion:**

1. **No free parameters beyond $\alpha$.** The entropy production criterion requires an additional coupling $\beta$; the structured record criterion requires an arbitrary persistence threshold $\tau_{\text{threshold}}$.
2. **Consistent with photon decoupling.** Photons don't decohere in vacuum. The CMB photon bath is already at thermal equilibrium. Only massive degrees of freedom under environmental monitoring contribute.
3. **Physically grounded.** Decoherence is measurement in the decoherence program. The framework claims duration requires measurement. The source of metering should be the rate of measurement events.
4. **Correct hierarchy:**

| Environment | $\gamma_D$ (s$^{-1}$) | $\rho_{\text{states}}$ (m$^{-3}$) | $J$ (m$^{-3}$ s$^{-1}$) |
|-------------|------------------------|--------------------------------------|---------------------------|
| Lab (condensed matter) | $10^{12}$ | $10^{28}$ | $10^{40}$ |
| Stellar core | $10^{15}$ | $10^{31}$ | $10^{46}$ |
| Molecular cloud | $10^{3}$ | $10^{10}$ | $10^{13}$ |
| Diffuse ISM | $10^{-5}$ | $10^{6}$ | $10^{1}$ |
| IGM | $10^{-15}$ | $1$ | $10^{-15}$ |
| Cosmic void (baryons) | $10^{-20}$ | $0.1$ | $10^{-21}$ |
| CMB photon field | $0$ | $0$ (massless) | $0$ |

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

**3. Executed Casimir benchmark null:** The repository now carries an explicit laboratory workbench in [`laboratory_bounds.py`](./laboratory_bounds.py), with the canonical report in [`laboratory_reports/casimir-benchmark-report.json`](./laboratory_reports/casimir-benchmark-report.json). Using the same Casimir-scale benchmark amplitude as above, together with quadratic weak-coupling scaling, the strongest current pressure-null anchor in the report is the Decca et al. 2003 parallel-plate measurement at `0.2 um`, which gives

$$\alpha \leq 6.79 \times 10^{-3}.$$

The secondary Decca et al. 2007 relative-error anchor gives

$$\alpha \leq 1.09 \times 10^{-2}.$$

At the `0.2 um` reference separation used in the report, this benchmark model implies:

- $\alpha = 10^{-2}$ would require an anomalous pressure of `1.30 mPa`
- $\alpha = 10^{-3}$ would require an anomalous pressure of `13.0 uPa`

So the first hard laboratory cutoff now exists in the repository, and it is already stronger than the older internal WEC consistency cutoff.

**4. Dimensional analysis:** In SI units, $\alpha_{\text{SI}} \sim \ell_P^3 \sim 4 \times 10^{-105}$ m$^3$/bit. In a laboratory ($\mu \sim 10^{26}$ bits/m$^3$): $\alpha_{\text{SI}} \cdot \mu \sim 4 \times 10^{-79}$.

**Interpretation:** The hardware test supports the lattice-model claim that a $1/N^2$ potential can produce the observed bound-state spectral structure. The Casimir workbench adds the first executed measurement-anchored laboratory null bound in the repository. But the lattice encodes $\alpha$ in circuit parameters, and the current Casimir bound still inherits the repository's benchmark anomaly model rather than a completed first-principles laboratory derivation. A direct physical extraction of spacetime $\alpha$ remains open.

---

## Metering Field Equation

### Action and Dynamics

The metering field $\mu(x, t)$ has its own dynamics:

$$S_\mu = \int \left[ \frac{(\partial_t \mu)^2}{2N} - \frac{N(\partial_x \mu)^2}{2} - N \cdot V(\mu) + N \cdot J \cdot \mu \right] dx \, dt$$

The Euler-Lagrange equation:

$$\partial_t^2 \mu = N^2 \partial_x^2 \mu + \frac{N'}{N}(\partial_t \mu)^2 + N N'(\partial_x \mu)^2 - N^2 V'(\mu) + N^2 J$$

**Propagation speed of $\mu$:** $c_\mu = N(\mu)$. In metered regions ($\mu \gg 1$), $c_\mu \to 1$. In unmetered regions ($\mu \to 0$), $c_\mu \to 0$. The metering density obeys the same causal structure as matter. Horizons are self-reinforcing: once $N = 0$, no $\mu$ dynamics can penetrate.

### Steady State

In the static limit with the free massive potential $V(\mu) = m^2 \mu^2 / 2$:

$$-\nabla^2 \mu + m^2 \mu = J(x)$$

This is the screened Poisson equation. The screening length $\ell_s = 1/m$ determines how far metering influence extends beyond the source.

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

Massive neutrinos ($m_\nu \sim 0.05$--$0.5$ eV) couple to the metering metric. Photons decouple and propagate at $c$. The excess delay for a massive particle traversing a void:

$$\Delta t = \int_{\text{void}} \left(\frac{1}{N(x)} - 1\right) \frac{dx}{c}$$

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

**Observational target:** Galaxy cluster lensing residuals. The archival program tests whether the sign-flip-like morphology survives model-family variation and whether any decoherence-relevant proxy carries conditional signal beyond stellar mass and local convergence. That is a morphology-and-confounder program, not yet a completed derivation of a direct photon-sector effect.

The current archival program works on the two public HFF cluster fields (`Abell 370`, `RXC J2248 / Abell S1063`) using the official HFF property catalogs, BUFFALO v2 photometric proxy catalogs, public MUSE Lensing Clusters spectroscopic catalogs, a hardened `radial_median_bandpass` residual construction, deeper external spectroscopy for membership hardening, public `MUSE-DEEP` core-cube extractions, and 22 independent frontier-model `kappa` maps. The sign-flip-like morphology is reproducible across model families, so the observational question is no longer whether the morphology exists in the archive; it is which explanatory axis survives the controls and replicates across clusters.

At present, `local_kappa` is the strongest and most stable predictor. Across the 22 model-family runs it is positive in `21/22`, significant in `19/22`, and has median partial correlation `0.3604`. Stellar mass remains a real axis as well: positive in `19/22`, significant in `16/22`, with median partial correlation `0.3740`. In the primary two-cluster pass, both combine at `p \approx 5.35 \times 10^{-5}`.

The BUFFALO-based continuum proxy does not survive the hardened residual. Across the 22 model-family runs it is positive in `8/22`, significant in `1/22`, with median partial correlation `-0.0316`.

The stronger spectroscopic pass is more informative, but it does not yet produce replication. MUSE-derived line-count and optical-complexity proxies are positive in `Abell 370` and negative in `RXC J2248 / Abell S1063` in the primary two-cluster pass. Across the 22 model-family runs, line count is positive in `14/22` and significant in `7/22`, with median partial correlation `0.1002`; optical complexity is positive in `14/22` and significant in `2/22`, with median partial correlation `0.1040`; emission strength is positive in `15/22` and significant in `2/22`, with median partial correlation `0.0786`. Those are real non-null structures, but not same-sign replication across the two clusters.

The first direct strong-line metallicity pass has now also been carried through using public MUSE line catalogs. In `Abell 370`, `N2` and Marino et al. (2013) oxygen-abundance estimators are mildly positive but not significant once mass, redshift, magnification, radius, and `local_kappa` are controlled. In `RXC J2248 / Abell S1063`, the same estimators are effectively sample-limited under the current public coverage. Across the 22 model-family runs, `muse_n2` is positive in `6/22` and significant in `0/22` with median partial correlation `0.0304`; the Marino `O3N2` abundance estimator is positive in `7/22` and significant in `0/22` with median partial correlation `0.0672`; the Marino `N2` abundance estimator is positive in `4/22` and significant in `0/22` with median partial correlation `-0.0480`.

The public matched-depth core-cube route has now been tested as well. Secure-redshift aperture extractions from the public `MUSE-DEEP` core cubes do not rescue the metallicity lane: the cubes cut off near `8850-8900 A`, so `Halpha/[N II]` fall out of band or effectively out of band at the cluster redshifts, and the remaining blue-line support is too thin to carry a two-cluster abundance test. In the current deep-cube extraction cache, `Abell 370` yields only `2` usable `R23` objects and `9` usable `O32` objects, while `RXC J2248 / Abell S1063` yields `0` usable `R23` objects and `1` usable `O32` object; `N2`, `O3N2`, and Marino-style abundance proxies remain unsupported in both clusters. That closes the observational proxy lane on current public data. Two later observational pursuits remain: redder symmetric spectroscopy across both clusters with usable `Halpha`, `[N II]`, `Hbeta`, and `[O III]`, and direct per-galaxy abundance tables mapped one-to-one into the HFF member catalogs. The active center of effort now shifts to the 3+1D lensing derivation.

### GW170817 Compatibility

GW170817/GRB 170817A constrains $|c_{\text{GW}} - c_{\text{photon}}|/c < 10^{-15}$.

Since the metering coupling enters through the mass term, all massless fields (photons, gravitons, gluons) decouple. Gravitational waves propagate at $c$ everywhere, regardless of $\mu$. Within this mass-term-only coupling model, the GW170817 constraint is satisfied for any $(m, \alpha)$, with no parameter tuning required.

---

## Observational Proposals

### A. Zeno Spectral Test (Completed)

A quantum system under spatially varying measurement strength exhibits discrete bound-state spectral structure matching the $1/N^2$ potential, rather than generic Zeno freezing. Realizable with current ion trap or superconducting qubit technology. **Three IBM Quantum hardware runs completed** (Section: [IBM Quantum Hardware](#ibm-quantum-hardware)).

### B. Anomalous Casimir Effect

The WEC violation predicts an additional negative-energy contribution in regions of high metering density. That lane is no longer just a proposal. A first executable laboratory workbench now exists in [`laboratory_bounds.py`](./laboratory_bounds.py), with the canonical report in [`laboratory_reports/casimir-benchmark-report.json`](./laboratory_reports/casimir-benchmark-report.json).

Under the repository's current benchmark model:

- the strongest executed Casimir pressure null gives $\alpha \leq 6.79 \times 10^{-3}$
- a signal at $\alpha = 10^{-2}$ would require an active-vs-inert pressure difference of about `1.30 mPa` at `0.2 um`
- a signal at $\alpha = 10^{-3}$ would require about `13.0 uPa` at the same separation

So this remains the shortest current path to a hard physical bound, but it is no longer vague. The next laboratory move is a purpose-built active-vs-inert Casimir differential with sub-`10 uPa` reach if the target is to penetrate below the `10^-3` regime in the present benchmark family.

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

The Einstein-backreaction branch is in better shape. In the exterior it is a GR-derived branch and relaxes toward Schwarzschild with $\gamma = 1$. That makes it the current mainline candidate. The new theory result is now six-stage. In spherical cluster-style maps, the branch remains convergence-dominated and does not produce a resolved raw sign-flip morphology. In toy multi-component composite maps with the same `radial_median_bandpass` style residual used observationally, the branch does produce robust anisotropic sign-changing residual structure while keeping the raw convergence map nonnegative. In a first real-member calibration built from the hardened HFF member tables, that same residual behavior appears in both `Abell 370` and `RXC J2248 / Abell S1063`. In the next direct archive-vs-theory step, a broad smooth envelope built from the same member geometry improves the structural comparison on the full 22-model ensemble: median resolved-sign agreement reaches `\approx 0.497`, median negative-residual overlap reaches `\approx 0.258`, median member-score Spearman reaches `\approx 0.165`, and median member pass-agreement reaches `\approx 0.658`. A subsequent rigid-alignment audit improves those structural and member-level comparisons further, especially the member-score ranking, but still leaves the field-level Pearson correlation near zero. So the immediate theoretical task has changed again. It is no longer to invent a photon channel, and no longer merely to ask whether multi-component geometry matters. It is to calibrate the viable Einstein branch beyond the first member-plus-envelope stage until the direct residual-map comparison is quantitatively competitive in the actual field, not only in sign structure and member-level rankings.

### The Coupling Constant

$\alpha$ is the central unknown. The hardware fit ($\alpha \approx 1.16$) applies to the lattice model, not to spacetime. The internal WEC benchmark ($\alpha < 0.133$) is a consistency cutoff, not a laboratory result. The repository now also contains a first executed measurement-anchored Casimir null bound, with strongest current benchmark value $\alpha \leq 6.79 \times 10^{-3}$ under the present anomaly model. The dimensional estimate ($\alpha_{\text{SI}} \sim \ell_P^3$) still suggests the physical coupling may be extraordinarily small. All cosmological predictions depend on $\alpha$, and a direct first-principles physical extraction of $\alpha$ remains open.

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

The mass-term coupling means $m_{\text{eff}} = m_0 / N(x)$. Heavier particles have larger absolute mass divergence in voids. The species hierarchy for void propagation is: heavy particles freeze first, light particles last. Neutrinos (lightest massive) are least affected. This may widen the detectable parameter window if the critical $N$ threshold varies by species.

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
| 3+1D lensing observable | Flat control, Einstein-backreaction branch, and explicit photon-sector branch are computable; the direct photon branch is heavily constrained; the Einstein branch is the current mainline candidate; spherical cluster-style maps remain convergence-dominated; multi-component, HFF-member, and calibrated member-plus-envelope Einstein maps generate sign-changing residual morphology without raw negative convergence; direct archive comparison now shows improved sign-structure and member-level agreement, and the rigid-alignment audit shows part of the remaining gap is registration-sensitive, but pixelwise field matching remains weak |
| Decoherence rate source criterion | Proposed and computationally explored; not observationally established |
| Void metering + cosmological predictions | Numerical parameter study complete; observational status open |
| Archival lensing program (2 clusters, 22 model families, BUFFALO + MUSE + deeper external spectroscopy + public MUSE-DEEP core cubes) | Morphology replicated under the hardened residual; local convergence and mass remain stable; the public-data proxy lane is exhausted without a replicated abundance-style signal; active focus shifts to the 3+1D lensing derivation |
| Coq formalization | Initial abstract theorem development complete; monolithic file compiles |
| Physical $\alpha$ laboratory bound | First Casimir benchmark null executed; strongest current benchmark cutoff $\alpha \leq 6.79 \times 10^{-3}$; direct first-principles measurement remains open |

---

## License

MIT
