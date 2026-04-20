# Quantum Metrology Benchmarks under Noise

This document summarizes the formulas used for calculating the Quantum Cramér-Rao Bound (QCRB) for both Holland-Burnett (HB) states and N00N states under photon-loss noise.

## 1. Noise Model: Independent Photon Loss
We assume each of the $d$ signal modes undergoes independent photon loss with transmission probability $\eta$ (where $\eta=1$ is ideal). The reference mode (mode 0) is assumed to be lossless.

$$E_k = \sum_{n=k}^{\infty} \sqrt{\binom{n}{k}} \eta^{(n-k)/2} (1-\eta)^{k/2} |n-k\rangle\langle n|$$

> [!NOTE]
> **Transition from Project A**: While Project A used direct **Stinespring dilation** (modeling environmental modes explicitly), this project uses **Kraus operators**. The two are mathematically identical, but the Kraus representation is significantly more memory-efficient as it avoids doubling the number of modes. This optimization is what makes computing the QCRB for $N=12$ and $N=18$ possible on a standard machine.

---

## 2. Holland-Burnett (HB) States (Numerical)
For HB states, the QCRB is computed numerically as there is no simple closed-form analytical expression for the noisy QFIM.

1.  **State Preparation**: The pure HB state $|\psi_{HB}\rangle$ is initialized.
2.  **Noise Application**: The density matrix $\rho_{out}$ is computed by applying the Kraus channel to all signal modes.
3.  **Quantum Fisher Information Matrix (QFIM)**: The elements $F_{ab}$ are calculated using the spectral decomposition of $\rho_{out} = \sum_n \lambda_n |n\rangle\langle n|$:
    $$F_{ab} = \sum_{n,m} \frac{2}{\lambda_n + \lambda_m} \text{Re}\left( \langle n | \partial_a \rho | m \rangle \langle m | \partial_b \rho | n \rangle \right)$$
    where $\partial_a \rho = -i[\hat{n}_a, \rho]$ and $\hat{n}_a$ is the number operator for mode $a$.
4.  **Total Variance**: The bound is the trace of the inverse QFIM:
    $$V_{HB} = \text{Tr}(F^{-1})$$

---

## 3. N00N State Benchmark (Proposed Analytical)
To compare HB states with N00N states under the same noise, we use $d$ independent N00N states, each with $N' = N/d$ photons, to estimate the $d$ phases.

### Single Phase N00N QFI (Asymmetric Loss)
For a single N00N state $\frac{1}{\sqrt{2}}(|N'\rangle_S |0\rangle_R + |0\rangle_S |N'\rangle_R)$ where mode $S$ (signal) has loss $\eta$ and mode $R$ (reference) is lossless, the QFI is:
$$F_{single}(\eta, N') = \frac{2 (N')^2 \eta^{N'}}{1 + \eta^{N'}}$$

### Total Variance for $d$ Phases
Summing the variances $V_j = 1/F_{single}$ for all $d$ phases, and substituting $N' = N/d$:
$$V_{N00N}(\eta, N, d) = \frac{d^3 (1 + \eta^{N/d})}{2 N^2 \eta^{N/d}}$$

**Verification**:
- If $\eta = 1$ (Noiseless), $V_{N00N} = \frac{d^3 (2)}{2 N^2 (1)} = \frac{d^3}{N^2}$. (Matches the ideal benchmark).
- If $\eta < 1$, the variance increases exponentially with $N$ due to the $\eta^{-N/d}$ term, reflecting the extreme sensitivity of N00N states to loss.

---

## 4. Summary Table

| State | Ideal Variance ($\eta=1$) | Noisy Variance ($\eta < 1$) |
| :--- | :--- | :--- |
| **Multipixel** | $\frac{d(1+\sqrt{d})^2}{4N^2}$ | Numerical (to be implemented if needed) |
| **N00N** | $\frac{d^3}{N^2}$ | $\frac{d^3 (1 + \eta^{N/d})}{2 N^2 \eta^{N/d}}$ |
| **HB** | Numerical | Numerical |
