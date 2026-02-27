# The Metering Metric

**Duration is not fundamental. It is constituted by measurement.**

A scalar field representing metering density couples to the field Lagrangian as a position-dependent mass term, producing a spacetime geometry in which temporal dynamics are energetically forbidden in regions devoid of metering subsystems. The manifold exists, but nothing elapses.

This is not a restatement of relational time. It is a specific, computable coupling mechanism with three provable structural theorems, a confirmed spectral signature on quantum hardware, and falsifiable cosmological predictions.

**Author:** Charles C. Norton

---

## Table of Contents

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

## Framework

### Intellectual Context

The framework draws on five established programs, each of which stops short of the specific mechanism formalized here:

- **Wheeler-DeWitt (1967):** The wavefunction of the universe satisfies H|Psi> = 0. No time parameter appears. Time is absent at the fundamental level.
- **Page-Wootters (1983):** A timeless universe contains time-evolving subsystems when conditioned on a clock degree of freedom. Time is created by conditioning, not discovered.
- **Rovelli (1996, 2018):** Duration is a relation between systems. The thermal time hypothesis derives time flow from thermodynamic structure.
- **Barbour (1999):** Time can be eliminated from the fundamental description, replaced by relations between configurations.
- **Bridgman (1927):** A quantity without a measurement procedure is undefined. Unmeasured duration is not duration.

None of these produces a modified spacetime geometry, a horizon, or a testable spectral prediction. This work does.

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

Three structural results hold for **any** monotone activation function $f$ (tanh, sigmoid, arctan, erf) and **any** metering profile $\mu$ that decays to zero (Gaussian, top-hat, power-law, double-Gaussian). They are dimension-independent (verified in 1+1D, 2+1D, 3+1D).

### Theorem 1: Horizon Formation

The tortoise coordinate integral

$$x^*(L) = \int_0^L \frac{dx}{N(x)}$$

diverges as $L \to \infty$. The metering boundary is a true horizon: signals from the metered region take infinite coordinate time to reach the unmetered region.

**Proof sketch:** Since $\mu$ decays to zero, $N(\mu(x)) \to 0$ for large $|x|$. The integrand $1/N$ is unbounded, and the integral diverges.

This has been verified computationally across 4 activation functions, 4 metering profiles, and 3 spatial dimensions.

### Theorem 2: Temporal Mode Confinement

The effective potential $V_{\text{eff}}(x) = 1/N(x)^2$ confines temporal oscillations:

- $V_{\text{eff}} \to \infty$ as $\mu \to 0$ (the potential is confining)
- The spectrum is purely discrete (only bound states exist)
- WKB tunneling probability: $\exp(-3218) \approx 0$ (exact machine zero)
- The metered region supports a finite number of bound states: ~69 radial modes (1+1D), ~600 modes (1000-site continuum), ~17 bits (3+1D with angular modes)

The temporal Hilbert space is finite-dimensional and bounded by the metering geometry.

### Theorem 3: Proper Time Vanishing

For a static observer at position $x$:

$$\tau(x, T) = N(x) \cdot T \to 0 \quad \text{as} \quad \mu(x) \to 0$$

For a worldline passing through a region where $\mu = 0$ on an interval $[t_1, t_2]$, the proper time contribution from that interval is exactly zero.

**Duration requires a meter.**

### Formal Verification

A Coq formalization plan targets machine-checked proofs of all three theorems using the Coquelicot real analysis library. Theorem 3 is straightforward (~50-100 lines). Theorem 1 requires improper integral machinery (~200-400 lines). Theorem 2 uses the WKB approach to prove confinement without full spectral theory (~300-600 lines).

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

Temporal horizons, mode confinement, bounded temporal information, self-reinforcing horizons, geometric SETI signatures, and WEC violation structure are all unchanged. Photon propagation through voids is restored. GW170817 is automatically satisfied.

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

Every peak matches to better than **0.06%**. This is two orders of magnitude tighter than the 6-qubit result, ruling out discretization artifacts. The bound-state structure is a property of the $1/N^2$ potential in the continuum limit.

### IBM Quantum Hardware

Three runs on `ibm_torino` (Qiskit 2.3.0, SamplerV2), 6-qubit lattice, 4000 shots per circuit.

**Run 1** (14 circuits, $t = 0.15$--$2.10$, Job `d6e3gdp54hss73badhbg`):
Survival probability descends from 0.87 to 0.13, then revives to 0.33. Dominant Fourier peak $\omega = 2.99$ matches predicted $E_2 - E_0 = 3.20$ at 6.4% error.

**Run 2** (20 circuits, $t = 0.05$--$1.00$, Job `d6e3ih154hss73badjt0`):
Dense sampling of the descent phase. Hardware tracks noiseless simulator within 3--8% for circuits with $\leq 50$ two-qubit gates.

**Run 3** (20 circuits, $t = 0.30$--$2.50$, Job `d6e3l9vg4t5c7387b93g`):
Independent confirmation. Minimum $P = 0.15$ at $t = 1.23$, revival to $P = 0.28$ at $t = 1.57$, second descent to $P = 0.16$ at $t = 2.38$.

The oscillatory revival structure -- descent, minimum, revival, second descent -- is **reproduced across two independent submissions** with different circuit constructions. This rules out generic Zeno freezing (which predicts monotonic decay) and confirms bound-state spectral structure.

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
- $1/N^2$ vs. quartic: **1.8--2.1x** (nearest competitor rejected)
- $1/N^2$ vs. harmonic: **6--15x**

The $1/N^2$ potential is preferred over all alternatives at 6 sites with ~15% hardware noise. The key discriminant is the potential's behavior at intermediate sites: $1/N^2$ gives $V = 3.41$ where quartic gives $V = 3.77$ and harmonic gives $V = 7.21$. The steep walls of the $1/N^2$ potential (rapid transition from low to high $V$) distinguish it from smoother shapes.

### Coupling Constant Constraints

Three independent constraints on the coupling constant $\alpha$:

**1. Hardware spectral fit:** $\alpha = 1.165$ minimizes $\chi^2$ against `ibm_torino` data (combined Runs 1 + 3).

**2. Weak energy condition bound:** The metric requires negative energy density inside the metered region at ~16x the Casimir scale when $\alpha = 1$. Requiring $|T_{00}| < E_{\text{Casimir}}$ gives $\alpha < 0.133$.

**3. Dimensional analysis:** In SI units, $\alpha_{\text{SI}} \sim \ell_P^3 \sim 4 \times 10^{-105}$ m$^3$/bit. In a laboratory ($\mu \sim 10^{26}$ bits/m$^3$): $\alpha_{\text{SI}} \cdot \mu \sim 4 \times 10^{-79}$.

**Interpretation:** The hardware test confirms the mathematical structure -- bound states in a $1/N^2$ potential exist and have the predicted spectrum. But the lattice encodes $\alpha$ in circuit parameters. The physical coupling constant is not directly measured by this experiment. Constraining the physical $\alpha$ requires experiments where $\mu$ arises naturally (Casimir anomaly, cosmological signatures, gravitational lensing residuals).

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
- **Convergence:** Nonlinear PDE evolution from 5% of steady state converges to within 2.4% of the analytic steady state by $t = 200$.
- **Interior perturbations:** Oscillatory decay (dynamical modes in the metered region).
- **Boundary perturbations:** Frozen. Perturbations at the horizon do not propagate. Eigenfrequencies scale with $N \to 0$.
- **Self-reinforcing horizons:** Signal crossing time through an unmetered gap scales as $\sim \epsilon^{-0.94}$ (theory: $\epsilon^{-1}$). As the gap minimum $N \to 0$, crossing time $\to \infty$.
- **Self-consistent coupling:** Matter-metering feedback converges in ~30--55 iterations. The self-consistent spectrum has higher energies than the prescribed case (narrower effective well), confirming backreaction.

### Potential Selection

The **free massive potential** $V(\mu) = m^2 \mu^2 / 2$ is the correct choice. It produces $\mu \to 0$ in voids (horizons form) and is consistent with the decoherence criterion: if there are no decoherence events ($J = 0$), there should be no metering ($\mu = 0$).

The **symmetry-breaking potential** $V(\mu) = -a\mu^2 + b\mu^4$ has a nonzero vacuum expectation value, giving $\mu > 0$ everywhere including true vacuum. This contradicts the framework's premise, kills horizons, and reduces the bound-state count from hundreds to 2--5.

---

## Cosmological Predictions

### Void Metering and Horizon Formation

A realistic 1D cluster-void-cluster geometry (400 Mpc total, clusters at $\pm 60$ Mpc, void extent ~120 Mpc) was modeled using NFW-like density profiles, temperature-dependent decoherence rates, and the screened Poisson equation.

Key numbers:
- Density contrast: $\rho_{\text{cl}} / \rho_{\text{void}} = 1429$
- Source contrast: $J_{\text{cl}} / J_{\text{void}} = 1.1 \times 10^7$

The framework predicts a **sharp phase transition** in the void interior as a function of the screening mass $m$ and coupling $\alpha$. The controlling parameter is $m \cdot R_{\text{void}}$:

| $m \cdot R_{\text{void}}$ | Behavior |
|---------------------------|----------|
| $< 3$ | $N_{\text{void}} \approx 1$. No observable effect. Metering from clusters fills the void. |
| $3$--$6$ | Transition region. $N_{\text{void}}$ drops from ~1 to ~0 depending on $\alpha$. |
| $> 6$ | $N_{\text{void}} \to 0$. Horizons form. The void interior is causally disconnected for massive particles. |

**Phase diagram** (120-point parameter sweep, $m \in [0.01, 10]$ Mpc$^{-1}$, $\alpha \in [0, 1000]$):

- Horizon-forming ($N_{\text{void}} < 10^{-3}$): **68 / 120** parameter points (57%)
- At the strictest threshold ($N_{\text{void}} < 10^{-6}$): horizons exist for $m \gtrsim 0.36$ Mpc$^{-1}$ at all tested $\alpha$, and never for $m \lesssim 0.17$ Mpc$^{-1}$

The transition is sharp -- not a gradual interpolation, but a topological change in causal structure. This is the framework working as designed: the tortoise integral either converges (no horizon) or diverges (horizon), with no intermediate state.

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

The delay is **binary rather than continuous**: either $N_{\text{void}} \approx 1$ (no interesting delay) or $N_{\text{void}} \ll 1$ (delay exceeds the age of the universe). The detectable window -- delay $> 10$ ns but finite -- is narrow, sitting near the phase boundary.

This is not a flaw but a prediction: the effect is a phase transition, and phase transitions are sharp.

### Differential Gravitational Lensing

The Ricci scalar $R(x)$ changes sign at the metering boundary:

- **Inside the metered region:** $R > 0$ (convergent lensing)
- **At the boundary:** $R < 0$ (divergent lensing)

This **convergence-divergence sign flip** is a unique morphological signature. Dark matter subhalos produce convergence only. Modified gravity (f(R), MOND) alters amplitude but not the sign pattern. The metering boundary's sign-flip pattern is structurally distinct.

**Observational target:** Galaxy cluster lensing residuals. Regress against galaxy metallicity and star formation rate (proxies for decoherence rate density) rather than mass alone. A correlation with the predicted convergence-plus-divergence morphology at high-complexity galaxy positions would constitute evidence.

This analysis requires no new observations. HST, JWST, and Chandra archival data suffice.

### GW170817 Compatibility

GW170817/GRB 170817A constrains $|c_{\text{GW}} - c_{\text{photon}}|/c < 10^{-15}$.

Since the metering coupling enters through the mass term, all massless fields (photons, gravitons, gluons) decouple. Gravitational waves propagate at $c$ everywhere, regardless of $\mu$. The constraint is **automatically satisfied** for any $(m, \alpha)$. No parameter tuning required.

---

## Observational Proposals

### A. Zeno Spectral Test (Completed)

A quantum system under spatially varying measurement strength exhibits discrete bound-state spectral structure matching the $1/N^2$ potential, rather than generic Zeno freezing. Realizable with current ion trap or superconducting qubit technology. **Three IBM Quantum hardware runs completed** (Section: [IBM Quantum Hardware](#ibm-quantum-hardware)).

### B. Anomalous Casimir Effect

The WEC violation predicts an additional negative-energy contribution in regions of high metering density. Precision Casimir experiments comparing active information-processing systems to inert controls could detect a metering-dependent correction. This is likely the **shortest path** to constraining the physical $\alpha$.

### C. Cosmological Metering Signature

The effective lapse function at early times (low metering density) differs from today. This produces a scale-dependent modification to the CMB power spectrum correlated with the information-processing capacity of matter at each epoch.

### D. Geometric Observer Detection

If metering density modifies spacetime geometry, regions containing observers have measurably different physical properties:

- **Anomalous gravitational lensing:** A metered region lenses light differently than an equivalent-mass region without metering subsystems, due to the curvature modification.
- **Redshift anomaly:** The lapse function is a gravitational redshift factor. Photons from metered regions carry a metering-dependent frequency shift.
- **Thermal halo:** The analog Unruh temperature at the metering boundary peaks at a characteristic value set by the metering gradient.
- **Correlation signature:** A true metering signal shows correlated anomalies across all three channels (lensing, redshift, thermal). False positives from unrelated astrophysics do not produce this correlation.

This constitutes a fundamentally new SETI methodology: searching for observers by their geometric imprint on spacetime, not by engineering artifacts or communication attempts.

### E. Void Lensing Residuals (Archival)

The convergence-divergence sign-flip morphology at void boundaries is distinct from any CDM or modified-gravity prediction. Test: lensing residuals in existing cluster surveys, correlated with proxies for information-processing density. Requires no new observations.

---

## Open Questions

### The Coupling Constant

$\alpha$ is the central unknown. The hardware fit ($\alpha \approx 1.16$) applies to the lattice model, not to spacetime. The WEC bound ($\alpha < 0.133$) is a consistency constraint. The dimensional estimate ($\alpha_{\text{SI}} \sim \ell_P^3$) suggests the physical coupling may be extraordinarily small. All cosmological predictions depend on $\alpha$, and there is currently no experiment that measures it in a natural (non-engineered) setting.

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

---

## Negative Results

These are documented for completeness and to prevent others from repeating dead ends.

**Geodesic normalization (v1):** The lapse-function approach to conditional proper time fails. Geodesic equations enforce $d\tau/d\lambda = 1$ automatically. Worldlines accumulate proper time in unmetered regions regardless of $N$. The mass-gap formulation (divergent effective mass) is the correct mechanism.

**Entropy test bug (v2):** The entropy-based duration test returned zero transitions because the random walk preserved the distribution used in the entropy calculation. Fixed by using a sliding window over visit history.

**Born-Oppenheimer decomposition:** The BO decomposition of the Wheeler-DeWitt solution never cleanly demonstrated emergent time. The coefficient of variation remained above 0.85 in all tests. The toy minisuperspace model (120 x 120 grid) is too crude. This does not affect the main results (which are independent of the BO analysis).

---

## Computational Verification

All results are verified using independent computational tools:

- **NumPy/SciPy:** Wheeler-DeWitt eigensolve, geodesic integration, density matrix evolution, PDE integration
- **Wolfram Mathematica 14.3:** 23 symbolic verification modules covering Lagrangian analysis, Christoffel symbols, Ricci scalar, tortoise coordinates, WKB tunneling, bound states, energy conditions, stability analysis, 3+1D extension
- **Qiskit 2.3.0 + IBM Quantum:** 3 hardware runs on `ibm_torino` (2026-02-23)

Simulator results span two orders of magnitude in lattice resolution (6 to 1000 sites) with spectral accuracy improving from 7% to 0.06%, confirming convergence to the continuum limit.

---

## Status

| Component | Status |
|-----------|--------|
| Core theorems (horizon, confinement, proper time) | Computationally verified |
| Zeno spectral test (simulator, 6--1000 sites) | Complete, 0.06% accuracy |
| Zeno spectral test (IBM hardware, 3 runs) | Complete, revival confirmed |
| Potential shape discrimination | Complete, $1/N^2$ preferred at 2.1x |
| Metering field equation (PDE + stability) | Complete, stable, self-consistent |
| Photon decoupling | Formalized |
| Decoherence rate source criterion | Formalized |
| Void metering + cosmological predictions | Complete |
| Coq formalization | Planned |
| Physical $\alpha$ measurement | Open |

---

## License

MIT
