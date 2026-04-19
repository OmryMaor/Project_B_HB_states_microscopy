# Holland-Burnett (HB) States: Math & Physics Behind the Implementation

This document breaks down exactly how the Holland-Burnett (HB) states are mathematically formulated in the `initialize_HB_state(N, d)` Python function, why this specific choice represents a direct generalization of the classical HB setup to multiple dimensions, and the foundational physics references.

### 1. Conceptual Background
The original **Holland-Burnett (HB) state** was proposed in 1993 [1] for 2-mode (single phase) interferometry. Instead of trying to generate highly non-linear N00N states (which are notoriously fragile and difficult to construct physically), Holland and Burnett realized that you can achieve the identical scaling $1/N^2$ Heisenberg variance bound trace limit proxy by simply injecting two identical Fock states into the two input ports of a standard 50:50 beam splitter. 

In our multi-parameter imaging scenario ($d=2$ phases $\implies d+1=3$ modes total, evaluated for $N=3$ photons), we utilize the theoretical generalization of this concept for multi-port interferometry (SU($d+1$)).

### 2. The Multi-Mode Generalization
To construct the equivalent of a "50:50 beam splitter" for $m = d+1$ arbitrary spatial modes, we rely on a symmetric unitary rotation known as the **Quantum Fourier Transform (QFT)** or multi-port beam splitter. It mixes all inputs symmetrically.

In your code script, this takes place iteratively:

1. **Calculate the mode distribution**: 
    If you allocate $N=3$ photons across $m = 3$ modes (2 target phases + 1 reference), the generalized HB approach requires injecting an equal number of identical photons into each input port. 
    $n_{\text{per port}} = N / m = 1$. 

2. **Construct the Multi-Port Symmetric Beamsplitter**:
   We define a phase factor $\omega = e^{i \frac{2\pi}{m}}$. The creation operators corresponding to the mixed/output modes $b_k^\dagger$ are superpositions of the independent mode creation operators $a_j^\dagger$:
   $$ b_k^\dagger = \frac{1}{\sqrt{m}} \sum_{j=0}^{m-1} \omega^{j \cdot k} \hat{a}_j^\dagger $$

3. **Construct the Initial State**:
   To generate the final HB state, we iteratively apply $n$ photons to each of the generated mixed modes across the vacuum state. 
   $$ |\psi_{HB}\rangle = \mathcal{N} \prod_{k=0}^{m-1} (b_k^\dagger)^n |0, 0, \dots, 0\rangle $$
   (Where $\mathcal{N}$ is a normalization constant resolved explicitly by `state.unit()` in the script).

### 3. Connection to Our Implementation Constraints

When $N=3, d=2$, you have exactly $m=3$ modes. Thus, $n=1$ photon is injected per mode.
When calculating this mathematically through the script's tensor operations, this results in the exact linear superpositions constructed in your plots:

$$ |\psi_{HB}\rangle = 0.82|1,0,0\rangle - 0.41|0,1,0\rangle - 0.41|0,0,1\rangle $$

Because there are practically no phase-sensitive nonlinearities deployed to strictly force pure correlations during multi-port beam combining, you'll observe in your plots that its evaluated $\operatorname{Tr}(\mathcal{F}^{-1})$ limits approach, but never identically match the absolute Multipixel $0.3238$ bounds reached by explicit optimization formulations (which is normal and physically expected), but rather offer a robust, easier-to-generate experimental baseline state capable of tracking similarly.

### 4. References 
1. **Original Dual-Mode HB States:** Holland, M. J., & Burnett, K. (1993). *Interferometric detection of optical phase shifts at the Heisenberg limit.* Physical Review Letters, 71(9), 1355.
2. **Multi-Parameter Limits:** Humphreys, P. C., et al. (2013). *Quantum Enhanced Multiple Phase Estimation.* Physical Review Letters, 111, 070403. (Discussing that while unoptimized generalized setups approach standard metrics stably, only optimized correlated arrays fully trace absolute bound convergences).
