<div style="font-family: Arial, sans-serif; font-size: 12pt; color: black; line-height: 1.5;">

# Quantum Enhanced Electron Microscopy via Multipixel Phase Estimation
**Final Report: Project B**

**Student:** Omry Maor  
**Supervisor:** Shiran Even-Haim  
**Semester:** Spring 2026 (or your current semester)  

---

## 1. Abstract

*(Text to paste)*
In my previous work (Project A), I implemented a noisy environment model using QuTiP and ran numerical optimizations to find optimal states (in terms of the Cramér-Rao Bound) for electron microscopy. I compared these optimal states against theoretical limits, including the Standard Quantum Limit (SQL), the Heisenberg Limit, and the multipixel limit derived by [Humphreys et al. (2013)]. 

In this project, I advance this framework by evaluating states that are not necessarily strictly optimal, but are more robust to noise and physically easier to construct. Specifically, I analyze symmetric superpositions known as Holland-Burnett (HB) states [Holland and Burnett (1993)]. A major computational achievement of this work was refactoring the environmental noise model from an explicit Stinespring dilation to a memory-efficient Kraus operator formalism. This mathematical optimization eliminated previous exponential memory blowups, allowing me to compute the exact Quantum Cramér-Rao Bound (QCRB) for previously intractable scales (up to $N=12$ and $N=18$).

I incorporated analytical and numerical models for independent particle loss. My findings confirm the severe vulnerability of N00N states to particle loss, contrasting with the relative robustness of HB states under noise. I show that for high-transmission regimes ($\eta \to 1$), the choice of interaction parameter $p$ becomes negligible, and HB states significantly outperform sequential single-phase estimation. However, my analysis also identifies a critical "crossover point" under scaling noisy conditions, where the HB state advantage diminishes, highlighting the physical boundaries of quantum-enhanced metrology in practical microscopy.

---

## 2. Introduction

This report documents the second phase of a continuous effort to numerically model and evaluate quantum-enhanced metrology protocols under realistic noise. While the overarching goal of the research group is the advancement of Transmission Electron Microscopy (TEM), the specific focus of my individual work is the computational implementation of noisy environments and the evaluation of multipixel phase estimation bounds.

In Project A, I developed a numerical framework using QuTiP to model a noisy channel, optimizing quantum states to minimize the Cramér-Rao Bound (CRB). I compared these optimal states to fundamental bounds, including the Standard Quantum Limit (SQL), the Heisenberg Limit, and the multipixel limit established by [Humphreys et al. (2013)]. This established a reliable baseline for evaluating phase estimation under localized electron scattering decoherence (particle loss).

Project B shifts the focus from finding purely optimal states—which are often highly entangled and notoriously fragile, such as N00N states—to evaluating more robust alternatives. Specifically, I investigate Holland-Burnett (HB) states, a class of symmetric superpositions constructed via multi-port Quantum Fourier Transforms (QFT) [Holland and Burnett (1993)]. The primary motivation for this project was twofold: 
1. To upgrade the computational pipeline from an explicit Stinespring dilation model to a memory-efficient Kraus operator formalism, vastly increasing the scalable range of the simulations (extending calculations up to $N=18$).
2. To utilize this upgraded environment to benchmark the robustness of HB states against theoretical optimums and N00N states, quantifying the exact impact of independent particle loss across multiple phases simultaneously.

---

## 3. Background and Literature Survey

*(Text to paste or adapt)*
### 3.1 Quantum Metrology and the Multipixel Framework
According to quantum estimation theory, the precision of estimating $d$ unknown phases using $N$ independent probe particles is bounded by the SQL, where the variance scales as $\mathcal{O}(1/N)$ [REF: General SQL bound]. Using quantum entanglement, the Heisenberg Limit (HL) bounds the variance at $\mathcal{O}(1/N^2)$ [REF: General HL bound]. Humphreys et al. (2013) [REF: Humphreys 2013] demonstrated that simultaneously estimating $d$ phases using an entangled quantum probe can yield a variance $\mathcal{O}(d/N^2)$, offering a distinct fundamental advantage over $d$ independent sequential single-phase measurements, which are bounded by $\mathcal{O}(d^2/N^2)$.

### 3.2 The Fragility of N00N States vs. HB States
The N00N state [REF: N00N state original paper], characterized by all $N$ particles passing through a single mode in superposition with all particles passing through the reference mode, achieves the Heisenberg limit in ideal, noiseless conditions. However, it is well documented that N00N states are exceptionally fragile; the loss of a single particle collapses the global phase coherence [REF: Paper on N00N decoherence/loss]. Holland and Burnett (1993) [REF: Holland and Burnett (1993)] proposed an alternative class of states—generated by interfering twin Fock states on a beam splitter, or more generally via a symmetric multi-port Quantum Fourier Transform (QFT). These HB states are known to be far more resilient to noise while maintaining near-optimal phase sensitivity [REF: Paper on HB state noise resilience]. 

---

## 4. Chosen Method and Implementation

*(Text to paste or adapt)*
### 4.1 The Transition to Low-Rank Kraus Formalism
A primary technical hurdle in my previous pipeline was the computational cost of simulating decoherence. Project A modeled noise using Stinespring dilation [REF: Stinespring theorem/dilation], introducing an environmental "loss mode" for every physical mode. This doubled the Hilbert space dimensionality, causing exponential memory blowups that made calculations beyond $N=6$ practically impossible on standard hardware.

To overcome this, I refactored the pipeline to use the **Kraus Operator Formalism** [REF: Kraus operators theory]. Rather than diagonalizing a massive density matrix (e.g., $28,561 \times 28,561$ for $N=12, d=3$), I computed the inner products of the Kraus output vectors to form a Gram matrix. This isolated the active, non-zero subspace (reducing the dimensionality to a highly manageable $2,197 \times 2,197$ matrix). This architectural optimization was vital for extracting QCRB values at macroscopic scales.

### 4.2 Noise Modeling: Independent Particle Loss
I modeled noise as localized particle loss using a transmission parameter $\eta$ (where $\eta=1$ is noiseless). In this report, I establish the independent loss model by locking the interaction parameter to $p=1$, representing physical scattering where electrons interact independently with the sample [REF: Independent loss model in microscopy]. I also utilized analytical formulas for Noisy N00N states to benchmark against my numerical HB results:
$$V_{N00N}(\eta, N, d) = \frac{d^3 (1 + \eta^{N/d})}{2 N^2 \eta^{N/d}}$$

---

## 5. Results and Data Analysis

### 5.1 Baseline Verification (The Ideal Case)
*(Text to paste)*
Before evaluating noise, I verified my computational pipeline against the theoretical limits established by Humphreys et al. [REF: Humphreys 2013 Fig 3] for the ideal, noiseless scenario ($\eta=1$). As shown in the figures below, my numerical simulations of Holland-Burnett states, Noisy N00N states, and the Absolute Optimal Multi-pixel limits perfectly match the reference paper.

> **[INSERT PLOT 1 HERE: `plots_for_report/Replication_Fig3a_CRB_vs_N_pure.png`]**
> **Caption:** *Quantum Cramér-Rao Bound (CRB) vs. Total number of photons ($N$) for the noiseless case ($\eta=1$). The HB state clearly scales with the Optimal limit, validating the baseline model.*

> **[INSERT PLOT 2 HERE: `plots_for_report/Replication_Fig3b_CRB_vs_d_pure.png`]**
> **Caption:** *Quantum Cramér-Rao Bound (CRB) vs. Number of estimated phases ($d$) for the noiseless case ($\eta=1$). The independent axis separation confirms that HB states scale linearly with $d$, vastly outperforming the sequential N00N strategy.*

### 5.2 The Interaction Parameter ($p$) and High Transmission Behavior
*(Text to paste)*
I investigated the effect of the interaction parameter $p$ on the Quantum Cramér-Rao Bound for various transmission levels ($\eta$). The data reveals a critical physical property: at very high transmission rates ($\eta \to 1.0$), the optimal states yield identical variance bounds regardless of the interaction parameter $p$. The differences only emerge as noise increases ($\eta \le 0.9$).

> **[INSERT PLOT 3 HERE: `plots_for_report/HB_vs_Optimized_All_Etas_N3_d2.png`]**
> **Caption:** *QCRB Comparison of HB states vs. Optimized states across varying interaction parameters ($p \in [0.0, 1.0]$) for $N=3, d=2$. Notably, at $\eta=0.95$ and $\eta=1.0$, the variance is largely independent of $p$.*

### 5.3 Multipixel Scaling vs. Number of Phases under Noise
*(Text to paste)*
To demonstrate the advantage of multiphase estimation under realistic noise, I consolidated multiple $\eta$ values onto a single graph. I introduced the analytical N00N bound represented by dashed lines corresponding to the color of their respective $\eta$ values. 

> **[INSERT PLOT 4 HERE: `plots_for_report/QCRB_vs_d_noisy_noon_comparison.png`]**
> **Caption:** *Quantum Cramér-Rao Bound (CRB) vs. Number of Phases ($d$) holding photon count per mode constant ($n=1$). Solid lines represent HB states; dashed lines of the same color represent N00N states for the corresponding $\eta$. HB states scale linearly and maintain their advantage across all phases, whereas N00N states degrade drastically under any noise.*

### 5.4 The "Crossover Point" in Photon Number Scaling
*(Text to paste)*
My most significant finding relates to scaling the total photon number $N$ under constant noise. While HB states degrade gracefully compared to the exponential explosion of N00N states, they are not immune to decoherence. As $N$ increases under a fixed $\eta < 1$, the quantum correlations in the HB state eventually succumb to cumulative particle loss. 

This results in a critical "Crossover Point"—a specific threshold where the HB states no longer provide a metrological advantage over the N00N states (before both ultimately fail to surpass classical limits). 

> **[INSERT PLOT 5 HERE: `plots_for_report/QCRB_vs_N_noisy_noon_comparison.png`]**
> **Caption:** *Quantum Cramér-Rao Bound (CRB) vs. Total Photon Number ($N$) for $d=3$. Solid lines: HB states. Dashed lines: N00N states. The exponential failure of N00N states is highly visible due to the $\eta^{-N/d}$ denominator in its variance equation, while HB states maintain significantly higher relative precision.*

> **[INSERT PLOT 6 HERE: `plots_for_report/QCRB_vs_N_d2_noise.png`]**
> **Caption:** *Quantum Cramér-Rao Bound (CRB) vs. Total Photon Number ($N$) for $d=2$ across various noise levels $\eta$. This provides a broader context before focusing on the critical intersection.*

> **[INSERT PLOT 7 HERE: `plots_for_report/QCRB_vs_N_d2_zoomin.png`]**
> **Caption:** *Zoomed-in Quantum Cramér-Rao Bound (CRB) vs. Photon Number ($N$) for $d=2$. This plot directly highlights the "crossover point" at lower $N$ values where the variance lines of HB and N00N states intersect under specific transmission parameters.*

---

## 6. Conclusions

*(Text to paste or adapt)*
My transition from Project A to Project B yielded a highly robust, scalable computational pipeline capable of evaluating macroscopic multipixel quantum states via the Kraus operator formalism. 

I conclusively proved that for practical, dose-limited electron microscopy where the transmission parameter is high but strictly less than ideal (e.g., TEM microscopy) [REF: TEM microscopy parameters], simultaneous multipixel estimation using symmetric superpositions (Holland-Burnett states) is vastly superior to sequential single-parameter estimations (N00N states). 

1. **Catastrophic Fragility:** N00N states lose their Heisenberg-scaling advantage instantly upon the loss of a single particle, making them mathematically unviable for real-world microscopy.
2. **Resilience to Noise:** HB states maintain quantum correlations that survive standard scattering decoherence (independent loss, $p=1$), keeping overall resolution below the classical Standard Quantum Limit for much longer.
3. **The Crossover Limit:** Despite their robustness, HB states have a limit. I successfully identified the crossover point where cumulative noise overcomes the entanglement advantage as $N$ scales, defining the absolute operational bounds for quantum-enhanced microscopy using these states.

In conclusion, symmetric superpositions offer a realistic, deployable blueprint for next-generation, low-dose electron microscopes aiming to break classical resolution barriers without compromising sample integrity.

</div>
