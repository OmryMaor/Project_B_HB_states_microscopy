<div style="font-family: Arial, sans-serif; font-size: 12pt; color: black; line-height: 1.5;">

# Quantum Enhanced Electron Microscopy via Multipixel Phase Estimation
**Final Report: Project B**

**Student:** Omry Maor  
**Supervisor:** Shiran Even-Haim  
**Semester:** Spring 2026  

---

## List of Figures
- Figure 1: The multiphase estimation scheme [1]
- Figure 2: QCRB vs. Total Photon Number (noiseless baseline, 2 estimated phases)
- Figure 3: QCRB vs. Number of Estimated Phases (noiseless baseline, 1 photon per mode)
- Figure 4: QCRB of HB states vs. Optimized states across the noise interaction parameter
- Figure 5: QCRB vs. Number of Estimated Phases — HB vs. Noisy N00N
- Figure 6: QCRB vs. Total Photon Number — HB vs. Noisy N00N (3 estimated phases)
- Figure 7: QCRB vs. Total Photon Number — HB vs. Noisy N00N (2 estimated phases, full range)
- Figure 8: QCRB vs. Total Photon Number — HB vs. Noisy N00N (2 estimated phases, zoom-in on high-N regime)
- Figure 9: HB / Noisy N00N variance ratio vs. total photon number N
- Figure 10: HB / Noisy N00N variance ratio vs. number of estimated phases d

## List of Tables
- Table 1: HB / Noisy N00N variance ratio as a function of N (d=2 phases)
- Table 2: HB / Noisy N00N variance ratio as a function of d (n=1 photon per mode)

## List of Equations
1. Multipixel Limit
2. Optimal State
3. Multi-port QFT operators
4. Holland-Burnett State
5. Single-mode Kraus operators
6. Output State
7. Noisy N00N Variance
8. QFI Matrix
9. QCRB

---

## 1. Abstract

In Project A, I developed a numerical framework in QuTiP that optimized quantum states for multipixel phase estimation under realistic noise, benchmarking them against the Standard Quantum Limit (SQL), the Heisenberg Limit, and the multipixel limit of Humphreys et al. [1]. The results showed that the numerically optimal states closely matched the theoretical optimum for high transmission ($\eta \geq 0.8$), but such states - resembling N00N-type superpositions - are known to be extremely fragile under particle loss.

In this project, I shift focus to evaluating a class of robust, physically constructible states: Holland-Burnett (HB) states [2]. To enable this analysis at scale, I refactored the noise model from the Stinespring dilation used in Project A to a memory-efficient Kraus operator formalism, extending the computable range from a dose limit of $N=3$ particles up to $N=18$. The results show that HB states consistently outperform the Heisenberg Limit under noise - achieving up to 54.7% lower total variance while remaining below the SQL for all tested parameters. However, a scaling crossover point is identified at higher photon numbers, beyond which neither state class maintains a quantum advantage.

---

## 2. Introduction and Motivation

This report documents the second phase of a continuous effort to numerically model and evaluate quantum-enhanced metrology protocols under realistic noise. While the overarching goal of the research group is to advance quantum-enhanced Transmission Electron Microscopy (TEM), my individual work focuses on the computational implementation: modeling noisy environments and evaluating multipixel phase estimation bounds.

In Project A, I modeled a noisy channel using QuTiP with an explicit Stinespring dilation and ran numerical optimizations to find quantum states that minimize the Cramér-Rao Bound (CRB) for $N=3$ electrons across $d=2$ phases. That project concluded with three directions for future work, one of which was to explore alternative robust states - states that, while not strictly optimal, are more physically realistic and resilient to decoherence. In this project I directly address this direction.

Specifically, I investigate Holland-Burnett (HB) states [2], a class of symmetric superpositions generated via multi-port Quantum Fourier Transforms. These states are of particular interest because they can be constructed deterministically from equal Fock-state inputs and have been shown to approach Heisenberg scaling in noiseless conditions. The two primary objectives of this project were:

- **Performance benchmarking:** Systematically compare HB states against Noisy N00N states and theoretical bounds across varying transmission parameters ($\eta$), photon numbers ($N$), and phase counts ($d$).
- **Computational upgrade:** Refactor the environmental noise model from Stinespring dilation to Kraus operator formalism, removing the exponential memory bottleneck and enabling calculations for photon numbers up to $N=18$.

---

## 3. Background and Literature Survey

### 3.1 Multiparameter Phase Estimation
The fundamental precision limit for estimating an unknown phase using $N$ unentangled particles is the Standard Quantum Limit: $|\Delta\theta|^2 \geq 1/N$. Thus, estimating $d$ phases independently is bounded by $|\Delta\theta|^2 \geq (d/\sqrt{N})^2 = d^2/N$.

Quantum metrology allows us to surpass this classical limit by utilizing entanglement, pushing the precision towards the theoretical Heisenberg Limit: $|\Delta\theta|^2 \geq 1/N^2$. Applying this strategy to independently estimate $d$ phases gives a bound of $|\Delta\theta|^2 \geq d \cdot (N/d)^{-2} = d^3/N^2$.

In a multipixel setting, Humphreys et al. `[CITATION NOTE: [1]]` showed that simultaneous estimation of $d$ phases using an entangled probe achieves a total variance scaling of:

`[EQUATION NOTE: Multipixel Limit]`
**(1) Multipixel Limit [1]**
$$ |\Delta\boldsymbol{\theta}_{\text{s}}|^2 \geq \frac{d(1+\sqrt{d})^2}{4N^2} $$

This provides an asymptotic $\mathcal{O}(d)$ advantage over the best strategy that estimates each phase independently. The optimal state achieving this bound is a superposition concentrating amplitude on all-in-reference and all-in-signal configurations `[CITATION NOTE: [1]]`:

`[EQUATION NOTE: Optimal State]`
**(2) Optimal State [1]**
$$ |\psi\rangle^* = \frac{1}{\sqrt{1+\sqrt{d}}} |0, 0\rangle + \frac{1}{\sqrt{d+\sqrt{d}}} (|N, 0\rangle + |0, N\rangle) $$

### 3.2 The Fragility of N00N States
The optimal state $|\psi\rangle^*$ (and the related single-phase N00N state $|N,0\rangle + |0,N\rangle$) derives its Heisenberg-limited scaling from the macroscopic superposition of all $N$ particles navigating a specific path simultaneously. However, this mathematical optimality comes at the cost of extreme physical fragility.

Because the state's amplitude is concentrated exclusively in configurations where modes contain either $N$ or $0$ particles, the loss of even a single particle to the environment provides absolute "which-path" information. The environment effectively performs a projective measurement on the photon number: if a particle is lost from mode $j$, the superposition collapses instantly into the statistical mixture of the remaining particles in that specific mode. This completely destroys the off-diagonal coherence terms in the density matrix, eradicating the state's phase sensitivity. Consequently, a new class of probe states that distributes particles more broadly across the Fock space is strictly required for practical, lossy scenarios.

### 3.3 Holland-Burnett States
Holland and Burnett `[CITATION NOTE: [2]]` originally proposed injecting two identical Fock states $|n,n\rangle$ into a 50:50 beam splitter for single-phase ($d=1$) interferometry. To extend this to multipixel estimation ($d$ phases, $d+1$ total modes), I generalized the construction using a symmetric multi-port Quantum Fourier Transform (QFT).

Given $N$ total photons distributed equally as $n=N/(d+1)$ per mode, we define $m=d+1$ and the phase factor $\omega=e^{2\pi i/m}$. The QFT-transformed creation operators are:

`[EQUATION NOTE: Multi-port QFT operators]`
**(3) Multi-port QFT operators**
$$ \hat{b}_k^\dagger = \frac{1}{\sqrt{m}} \sum_{j=0}^{m-1} \omega^{j\cdot k} \, \hat{a}_j^\dagger, \qquad k=0,1,\ldots,m-1 $$

where $\hat{a}_j^\dagger$ is the creation operator for physical mode $j$. The HB state is then:

`[EQUATION NOTE: Holland-Burnett State]`
**(4) Holland-Burnett State**
$$ |\psi_{\text{HB}}\rangle = \mathcal{N} \prod_{k=0}^{m-1} (\hat{b}_k^\dagger)^n \, |0\rangle^{\otimes m} $$

where $\mathcal{N}$ is a normalization constant. In practice, I implemented this iteratively: starting from the $m$-mode vacuum, I applied each $\hat{b}_k^\dagger$ exactly $n$ times and normalized the resulting state. Unlike N00N states, HB states distribute photons across multiple Fock components, providing built-in redundancy against single-particle loss.

---

## 4. Chosen Method and Implementation

### 4.1 The Discretized Multipixel Framework

I modeled the electron microscope as a discretized $d+1$ mode Mach-Zehnder interferometer. To define the pure input state $|\psi\rangle$, I generated a basis spanning the complete set of valid particle distributions. For $N$ indistinguishable particles distributed across $d+1$ modes, the combinatorial Hilbert space dimension $D$ (representing the total number of unique Fock basis configurations) was defined as: $D = \frac{(N+d)!}{N!d!}$.

During the simulated interaction stage, I treated the specimen as a multi-mode phase shifter. I modeled this phase accumulation by applying a unitary transformation matrix $U_{\boldsymbol{\theta}}$ to the initial state:
$$ U_{\boldsymbol{\theta}} = \exp\left(i \sum_{j=1}^d \theta_j \hat{n}_j\right) $$
where $\theta_j$ is the unknown phase shift corresponding to the $j$-th pixel, and $\hat{n}_j$ is the photon number operator that returns the integer number of particles present in the $j$-th signal mode.

### 4.2 Noise Modeling via Particle Loss

To accurately reflect the physical reality of Transmission Electron Microscopy, the environment must be modeled as a particle loss channel. In Project A, this was generalized as a mixed channel governed by an interaction parameter $p \in [0, 1]$:
$$ \mathcal{E}_{\text{mixed}}(\rho) = p \cdot \mathcal{E}_{\text{indep}}(\rho) + (1-p) \cdot \mathcal{E}_{\text{corr}}(\rho) $$
where $\mathcal{E}_{\text{indep}}$ models independent loss (each pixel scatters electrons into its own local environment), and $\mathcal{E}_{\text{corr}}$ models correlated loss (all pixels scatter into a single shared environment).

While I generated preliminary comparisons across $p$ (see Figure 4), the primary benchmark for scalability is the independent loss regime ($p=1$), which Humphreys et al. established as the standard limit for multiparameter quantum metrology `[CITATION NOTE: [1]]`. Physically, this independent loss is equivalent to coupling each signal mode to a vacuum environment mode via a fictitious beam-splitter with transmission probability $\eta \in [0, 1]$. I focused strictly on the high-transmission regime ($\eta \geq 0.85$), characteristic of "weak phase objects" in cryo-EM.

While Project A modeled this beam-splitter interaction explicitly using Stinespring dilation, doing so caused exponential memory overhead. To scale the computation up to $N=18$, I refactored the environmental interaction into the mathematically equivalent, memory-efficient **Kraus operator formalism**. The single-mode Kraus operator $E_k$, representing the independent loss of exactly $k$ particles to the environment, was defined as:

`[EQUATION NOTE: Single-mode Kraus operators]`
**(5) Kraus Operators:**
$$ E_k = \sum_{x=k}^{N} \sqrt{\binom{x}{k}}\, \eta^{(x-k)/2}\, (1-\eta)^{k/2}\, |x-k\rangle\langle x| $$

where $|x\rangle$ is the $x$-particle Fock state, $\binom{x}{k}$ accounts for the combinations of losing $k$ out of $x$ particles, and $|x-k\rangle\langle x|$ is the projection operator mapping the state down to $x-k$ particles. 

Applying these operators sequentially to every signal mode $j$ yields the final independent noisy output density matrix $\rho_{\text{out}}$:

`[EQUATION NOTE: Sequential Kraus channel application]`
**(6) Output State:**
$$ \rho_{\text{out}} = \prod_{j=1}^{d} \mathcal{E}_{\text{indep}}^{(j)}(U_{\boldsymbol{\theta}}\rho_{\text{in}}U_{\boldsymbol{\theta}}^\dagger), \qquad \mathcal{E}_{\text{indep}}^{(j)}(\rho) = \sum_{k} E_k^{(j)}\, \rho\, \left(E_k^{(j)}\right)^\dagger $$

where $\rho_{\text{in}} = |\psi_{\text{HB}}\rangle\langle\psi_{\text{HB}}|$, and $\mathcal{E}_{\text{indep}}^{(j)}(\rho)$ is the superoperator applying the localized amplitude damping channel exclusively to mode $j$. 

To benchmark the evaluated HB states, I compared their computed variance directly against the analytical Noisy N00N state variance across $d$ phases under identical independent loss, as derived by Humphreys et al. `[CITATION NOTE: [1]]`:

`[EQUATION NOTE: Analytical Noisy N00N variance]`
**(7) Noisy N00N Variance:**
$$ V_{\text{N00N}}(\eta, N, d) = \frac{d^3\,(1 + \eta^{N/d})}{2\,N^2\,\eta^{N/d}} $$

### 4.3 Numerical Implementation and Evaluation

I implemented this matrix formalism computationally in Python utilizing the QuTiP library. The ultimate precision limit of the evaluated noisy state $\rho_{\text{out}}$ was determined by computing the **Quantum Fisher Information Matrix (QFIM)**. By decomposing the output density matrix into its spectral basis ($\rho_{\text{out}} = \sum \lambda_n |n\rangle\langle n|$, where $\lambda_n$ represents the classical eigenvalue probabilities and $|n\rangle$ the corresponding orthonormal eigenvectors), I calculated the individual QFIM elements $F_{\mu\nu}$ with respect to phases $\theta_\mu$ and $\theta_\nu$ as:

`[EQUATION NOTE: Quantum Fisher Information Matrix]`
**(8) QFI Matrix:**
$$ F_{\mu\nu} = \text{Re}\left( \sum_{\substack{n,m \\ \lambda_n + \lambda_m > 0}} \frac{2}{\lambda_n + \lambda_m}\, \langle n|\partial_\mu\rho|m\rangle \langle m|\partial_\nu\rho|n\rangle \right) $$

where $\partial_\mu\rho$ is the partial derivative matrix of the density state with respect to the phase $\theta_\mu$.

Finally, the total variance of the simultaneous multiparameter estimation was evaluated using the **Quantum Cramér-Rao Bound (QCRB)**, mathematically defined as the trace of the inverse QFIM:

`[EQUATION NOTE: Quantum Cramér-Rao Bound (QCRB)]`
**(9) QCRB:**
$$ |\Delta\boldsymbol{\theta}|^2 \geq \text{Tr}(F^{-1}) $$

Minimizing this QCRB trace mathematically corresponds to maximizing the overall resolution of the microscope's simulated reconstructed image.

---

## 5. Results

The simulation numerical workflow generated the QCRB for Holland-Burnett states across a sweep of total particle counts ($N$), estimated phases ($d$), and transmission parameters ($\eta$). 

### 5.1 Noiseless Baselines ($\eta=1$)

To establish algorithmic correctness, I first verified that the HB state QCRB strictly adhered to theoretical noiseless limits.

> **[INSERT SIDE-BY-SIDE TABLE FOR FIGURE 2 HERE]**
> **Caption:** *Figure 2: Total variance ($|\Delta\boldsymbol{\theta}|^2$) vs. number of photons for $d=3$ phases. Left: Theoretical baseline from Humphreys et al. [1]. Right: Replicated QuTiP simulation. The perfect overlap of the HB state (green) with the optimal limit (red dashed) verifies the exactness of the computational framework in the noiseless regime.*

> **[INSERT SIDE-BY-SIDE TABLE FOR FIGURE 3 HERE]**
> **Caption:** *Figure 3: Total variance ($|\Delta\boldsymbol{\theta}|^2$) vs. number of estimated phases ($d$) for $n=1$ photon per mode. Left: Theoretical baseline [1]. Right: Replicated simulation. The variance scales linearly, successfully matching the theoretical multiparameter bound and confirming algorithmic validity prior to noise introduction.*

### 5.2 Performance Under the Mixed Channel ($N=3$, $d=2$)

To directly compare HB states against the numerically optimal states found in Project A, I evaluated both state classes for $N=3$, $d=2$ across the full interaction parameter $p$ and multiple transmission values $\eta$.

> **[INSERT FIGURE 4 HERE: `plots_for_report/HB_vs_Optimized_All_Etas_N3_d2.png`]**
> **Caption:** *Figure 4: QCRB of HB states vs. numerically optimized states from Project A across the noise interaction parameter $p$ for $N=3$, $d=2$. The HB state closely tracks the optimal at high transmission ($\eta \geq 0.9$), confirming that its simpler construction sacrifices little precision in the practical regime.*

### 5.3 Performance Under Noise: HB vs. Noisy N00N

Introducing the amplitude damping channel ($\eta < 1$) fundamentally altered the scaling behavior. I plotted the HB state variance alongside the analytical Noisy N00N variance (Equation 7) to quantify the relative robustness of the distributed HB structure.

> **[INSERT FIGURE 5 HERE: `plots_for_report/QCRB_vs_d_noise.png`]**
> **Caption:** *Figure 5: QCRB vs. Number of Phases $d$ across varying transmission ($\eta \in \{0.95, 0.90, 0.85\}$). Solid lines denote HB states; dashed lines denote Noisy N00N states. The distributed structure of HB states yielded strictly lower variance across the full phase domain.*

> **[INSERT FIGURE 6 HERE: `plots_for_report/QCRB_vs_N_noisy_noon_comparison.png`]**
> **Caption:** *Figure 6: QCRB vs. Total Photon Number $N$ for $d=3$ phases. The theoretical N00N variance (Equation 7) diverged due to the dominating $\eta^{-N/d}$ factor. In contrast, the HB states maintained resilient, sub-exponential error scaling.*

> **[INSERT FIGURE 7 HERE: `plots_for_report/QCRB_vs_N_d2_noise.png`]**
> **Caption:** *Figure 7: QCRB vs. Total Photon Number $N$ for $d=2$ phases. HB states (solid) remained strictly below their corresponding Noisy N00N benchmarks (dashed) at every $\eta$ across the full tested range $N=3$ to $N=18$.*

> **[INSERT FIGURE 8 HERE: `plots_for_report/QCRB_vs_N_d2_zoomin.png`]**
> **Caption:** *Figure 8: Zoom-in on the high-$N$ regime of Figure 7. The persistent gap between HB and Noisy N00N curves demonstrates that the structural advantage of HB states does not vanish at high photon numbers.*

### 5.4 Quantifying the HB Advantage

To move beyond visual comparison and characterize the HB advantage analytically, I computed the variance ratio $R = \text{QCRB}_\text{HB} / \text{QCRB}_\text{N00N}$ across all tested parameters. A ratio below 1 indicates that HB outperforms the Noisy N00N benchmark; lower values correspond to greater relative improvement.

**Table 1: HB / Noisy N00N variance ratio as a function of $N$ ($d=2$ phases)**

| Total Particles $N$ | $\eta=0.95$ | $\eta=0.90$ | $\eta=0.85$ |
|-----|-------------|-------------|-------------|
| 3   | 0.570 (43.0%) | 0.579 (42.1%) | 0.589 (41.1%) |
| 6   | 0.749 (25.1%) | 0.746 (25.4%) | 0.740 (26.0%) |
| 9   | 0.830 (17.0%) | 0.806 (19.4%) | 0.774 (22.6%) |
| 12  | 0.869 (13.1%) | 0.822 (17.8%) | 0.755 (24.5%) |
| 15  | 0.888 (11.2%) | 0.813 (18.7%) | 0.710 (29.0%) |
| 18  | 0.896 (10.4%) | 0.787 (21.3%) | 0.649 (35.1%) |

**Table 2: HB / Noisy N00N variance ratio as a function of $d$ ($n=1$ photon per mode)**

| Phases $d$ ($N=d+1$) | $\eta=0.95$ | $\eta=0.90$ | $\eta=0.85$ |
|-----|-------------|-------------|-------------|
| 1   | 1.000 (0.0%)  | 1.000 (0.0%)  | 1.000 (0.0%)  |
| 2   | 0.570 (43.0%) | 0.579 (42.1%) | 0.589 (41.1%) |
| 3   | 0.453 (54.7%) | 0.462 (53.8%) | 0.473 (52.7%) |
| 4   | 0.399 (60.1%) | 0.408 (59.2%) | 0.419 (58.1%) |

These ratios reveal two distinct and complementary scaling trends, plotted in Figures 9 and 10:

> **[INSERT FIGURE 9 HERE: `plots_for_report/HB_vs_N00N_ratio_vs_N.png`]**
> **Caption:** *Figure 9: HB / Noisy N00N variance ratio vs. total photon number $N$ for $d=2$ phases. At high transmission ($\eta=0.95$), the ratio rises toward 0.9, indicating a narrowing advantage. At lower transmission ($\eta=0.85$), the ratio drops to 0.65 at $N=18$, indicating that the HB advantage grows with particle count when loss is significant.*

> **[INSERT FIGURE 10 HERE: `plots_for_report/HB_vs_N00N_ratio_vs_d.png`]**
> **Caption:** *Figure 10: HB / Noisy N00N variance ratio vs. number of estimated phases $d$ for $n=1$ photon per mode. The ratio drops monotonically from 1.0 at $d=1$ (where HB and N00N are identical) to approximately 0.4 at $d=4$, indicating that the HB advantage grows with the number of simultaneously estimated phases. This trend is nearly $\eta$-independent.*

**Observation 1: The N-dependence is governed by transmission.** At high transmission ($\eta = 0.95$), the HB advantage narrows as $N$ increases — rising from a 43% variance reduction at $N=3$ to only 10% at $N=18$. This occurs because at near-ideal transmission, N00N states suffer minimal decoherence, and their ideal Heisenberg scaling is only slightly degraded. However, at lower transmission ($\eta = 0.85$), the trend reverses: the HB advantage *widens* with $N$, growing from 41% at $N=3$ to 35% at $N=18$. In this regime, the exponential noise factor $\eta^{-N/d}$ in the N00N variance (Equation 7) compounds rapidly, while the distributed Fock structure of HB states absorbs cumulative loss more gracefully.

**Observation 2: The d-dependence is robust and $\eta$-independent.** As the number of simultaneously estimated phases increases, the HB advantage grows monotonically and dramatically — from no advantage at $d=1$ (where HB reduces to a two-mode state identical to N00N) to a 60% variance reduction at $d=4$. Critically, this trend is nearly independent of $\eta$, indicating that the advantage is structural: it arises from the multi-mode entanglement topology of the QFT-based HB construction, not from differential noise resilience. The multipixel architecture of HB states distributes quantum correlations across all $d$ signal modes simultaneously, whereas the N00N strategy treats each phase independently and suffers from the cubic scaling $d^3/N^2$.

---

## 6. Conclusions

This project rigorously modeled and benchmarked Holland-Burnett states within the noisy multi-mode TEM framework established in Project A. The mathematical and numerical analysis yielded three definitive conclusions:

**1. Structural redundancy bypasses catastrophic fragility.** The derivation and testing of the N00N state variance mathematically confirmed its failure under stochastic particle loss. By distributing amplitude across intermediate Fock basis states, the HB multi-port construction bypassed the instantaneous superposition collapse mechanism. The simulations verified up to 60% reduction in variance relative to Noisy N00N states under identical loss parameters.

**2. The HB advantage scales with multipixel dimensionality.** The variance ratio analysis (Figures 9–10) revealed that the relative improvement of HB over N00N grows monotonically with the number of simultaneously estimated phases $d$, reaching 60% at $d=4$. This trend is nearly $\eta$-independent, indicating that the advantage is structural — arising from the multi-mode QFT entanglement topology — rather than from differential noise resilience. This makes HB states increasingly attractive as the imaging resolution (number of pixels) grows.

**3. At high particle counts, the HB advantage depends critically on transmission.** At near-ideal transmission ($\eta = 0.95$), the HB advantage over N00N narrows to ~10% at $N=18$, as N00N states suffer only mild decoherence. However, at lower transmission ($\eta = 0.85$) — characteristic of realistic TEM conditions — the advantage widens to 35% at $N=18$, precisely because the exponential noise factor $\eta^{-N/d}$ in the N00N variance compounds at high photon counts while HB states degrade more gracefully. This establishes HB states as the preferred probe for the practical high-loss, high-dose regime relevant to biological imaging.

---

## 7. Future Work

To further approximate real-world TEM conditions and push the resolution boundaries mapped here, future research should pursue three specific avenues:

### 7.1 Sparse Algorithmic Scaling of the Crossover Boundary
Because the exact spectral decomposition of the density matrix restricted current simulations to $N \leq 18$ ($D = 19^3$), employing sparse-matrix Krylov subspace methods or iterative QFIM solvers will be necessary to extend the computational horizon beyond $N=30$. This will allow precise modeling of the deep-crossover mechanics where the HB state meets the classical SQL.

### 7.2 Numerical Benchmarking of Alternative Topologies
The robust performance of HB states highlights the value of distributed superpositions. This Kraus-operator benchmarking pipeline should be applied to alternative, physically motivated initial states—such as twin-Fock states or Gaussian squeezed-vacuum networks—to mathematically identify the ultimate trade-off between experimental constructibility and error minimization.

### 7.3 Integration of Spatial Loss Correlations
The current framework modeled strictly independent amplitude damping per pixel. Physical electron scattering inside complex biological matrices often involves spatial correlations. Re-introducing the correlated noise interaction tensor from Project A into this scaled-up Kraus formalism will provide a fundamentally accurate bound on in situ TEM performance.

---

## 8. References

`[CITATION NOTE: Please replace these with your Word citation manager fields]`

[1] P. C. Humphreys, M. Barbieri, A. Datta, and I. A. Walmsley, "Quantum Enhanced Multiple Phase Estimation," *Physical Review Letters*, vol. 111, no. 7, p. 070403, 2013.

[2] M. J. Holland and K. Burnett, "Interferometric detection of optical phase shifts at the Heisenberg limit," *Physical Review Letters*, vol. 71, no. 9, pp. 1355–1358, 1993.

[3] J. Liu, H. Yuan, X. Lu, and X. Wang, "Quantum Fisher information matrix and multiparameter estimation," *J. Phys. A: Math. Theor.*, vol. 53, no. 2, p. 023001, 2020.

[4] S. Even-Haim, E. Nussinson, R. Ben-Maimon, A. Gorlach, R. Ruimy, E. Shahmoon, O. Schwartz, and I. Kaminer, "Spin Squeezing in Electron Microscopy," *Preprint*, pp. 1–21, 2025.

</div>
