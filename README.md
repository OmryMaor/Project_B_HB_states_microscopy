# Quantum Enhanced Electron Microscopy via Multipixel Phase Estimation

![Python](https://img.shields.io/badge/Python-3.10+-blue?style=for-the-badge&logo=python)
![QuTiP](https://img.shields.io/badge/QuTiP-Quantum_Toolbox-purple?style=for-the-badge)
![Jupyter](https://img.shields.io/badge/Jupyter-Notebook-orange?style=for-the-badge&logo=jupyter)

A computational physics project investigating simultaneous multipixel phase estimation in noisy environments, aimed at pushing the resolution limits of electron microscopes without increasing radiation dosage.

## 🔬 The Physics Problem
In electron microscopy, pushing the resolution limits often means exposing biological and beam-sensitive materials to higher electron radiation doses, causing irreversible sample damage. The goal of quantum metrology is to overcome the **Standard Quantum Limit (SQL)** by extracting more information per individual probing particle.

While my previous work (**Project A**) focused on developing a general numerical framework using QuTiP to optimize states under noise, **Project B** advances this framework. Rather than searching for theoretically optimal states that are heavily entangled and highly fragile (like N00N states), this repository focuses on benchmarking **Holland-Burnett (HB) states**. These states—generated via multi-port Quantum Fourier Transforms (QFT)—are significantly more robust to noise and physically easier to construct.

## 💻 Computational Achievements
A major hurdle in modeling quantum environments is the exponential memory explosion when simulating particle loss.

* **Legacy Approach (Stinespring Dilation):** Explicitly modeled the environment with "loss modes," doubling the physical modes and crashing computations beyond $N=6$ due to astronomical matrix dimensions.
* **Refactored Approach (Kraus Formalism):** I successfully refactored the pipeline to use a memory-efficient low-rank **Kraus Operator Formalism**. By extracting the density matrix derivatives via spectral decomposition and mapping them over Gram matrices, the exact Quantum Cramér-Rao Bound (QCRB) can now be computed for macroscopic scales up to **$N=18$**.

## 📊 Key Results
By modeling localized electron scattering as an amplitude damping channel (independent particle loss), I established the following bounds:

1. **Catastrophic Fragility of N00N States:** The N00N states lose their Heisenberg-scaling advantage instantly upon the loss of a single particle. 
2. **Robustness of HB States:** Holland-Burnett states maintain quantum correlations that survive standard scattering decoherence, degrading gracefully and keeping overall resolution below classical limits for significantly longer.
3. **The Crossover Point:** The simulation explicitly captures the intersection point where cumulative scaling noise eventually overcomes the multipixel entanglement advantage.

### High Transmission Robustness ($d=2$)
<img src="plots_output/plots_for_report/QCRB_vs_N_d2_noise.png" width="600"/>

> *At $\eta=1.0$, HB states scale identically to the optimal theoretical limit. Under noise ($\eta \le 0.95$), HB states remain far more robust than Noisy N00N counterparts.*

### The Crossover Point Detail
<img src="plots_output/plots_for_report/QCRB_vs_N_d2_zoomin.png" width="600"/>

> *Zoomed-in perspective highlighting the physical "crossover point" at lower $N$ values where the variance bounds of HB and N00N states intersect under specific transmission parameters.*

## 📁 Repository Structure
To keep this repository clean and readable for external AI analysis and review, all legacy scripts and data batches have been excluded. 
- `Project_B_Complete_Codebase.ipynb`: A chronologically ordered, deeply documented Jupyter Notebook containing the exact physics functions (Kraus loss, HB initialization, QFIM calculation) alongside a fast-execution plotting block.
- `Project_B_Comprehensive_Theory.md`: A highly readable, story-driven document integrating all mathematical derivations, equations, and theoretical context behind the codebase.
- `Project_B_Final_Report_Template.md`: The finalized academic report template detailing the formal findings.
- `plots_output/plots_for_report/`: The generated empirical plots.

## 🚀 Getting Started
The `Project_B_Complete_Codebase.ipynb` is designed to be plug-and-play. The physics engine is fully exposed for analysis, and the final plot generation blocks utilize pre-calculated matrices to instantly render the results without incurring a 45-minute eigen-decomposition overhead on local machines.
