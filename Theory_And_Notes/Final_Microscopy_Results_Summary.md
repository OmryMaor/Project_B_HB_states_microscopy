# Final Results Summary: Multipixel Phase Estimation under Kraus-modeled Decoherence

## 1. Core Computational Achievement: Low-Rank Kraus Formalism
We have successfully refactored the noise channel implementation, transitioning from an explicitly modeled Stinespring dilation (which suffered from exponential memory blowup) to a memory-efficient Kraus operator formalism. For our N=12 and d=3 evaluations, instead of eigendecomposing the full 28,561 × 28,561 density matrix, we isolated the output subspace by evaluating the Gram matrix of the composite Kraus output vectors. This reduced the active dimensional space to a 2,197 × 2,197 matrix, enabling the precise computation of the Quantum Cramér-Rao Bound (QCRB) for highly complex states that were previously computationally intractable.

## 2. Physical Insight: The Fragility of N00N States
Our results quantitatively validate the catastrophic vulnerability of N00N states to particle loss. We benchmarked the noisy N00N variance analytically using the derived formula:
$$V_{N00N}(\eta, N, d) = \frac{d^3 (1 + \eta^{N/d})}{2 N^2 \eta^{N/d}}$$
The presence of the transmission parameter $\eta$ in the exponential denominator proves that as the total photon number N or the number of phases d scales, the N00N state loses its Heisenberg-scaling advantage entirely. Even at a high transmission of $\eta=0.85$, the variance explodes because a single lost particle destroys the global phase coherence, reducing the state to a statistical mixture.

## 3. Physical Insight: The Robustness of Holland-Burnett (HB) States
The generated plots confirm that Holland-Burnett states, initialized via a symmetric Multi-port Quantum Fourier Transform (QFT), are exceptionally resilient to independent particle loss (where the mixture parameter is physically locked to $p=1$ to accurately model localized electron scattering in a microscope).

*   **CRB vs N Scaling**: As N scales up to 18, the HB states continue to demonstrate monotonic improvement. They maintain a tight variance boundary just above the absolute optimal multipixel limit, gracefully degrading rather than catastrophically collapsing like the N00N equivalents.
*   **CRB vs d Scaling**: When holding the per-mode particle count steady and increasing the number of estimated phases, the variance of the HB states scales almost linearly. They drastically outperform N00N states across all noisy regimes.

## 4. Strategic Conclusion for Quantum Microscopy
The data solidifies the core thesis of our research: for practical, dose-limited electron microscopy where the transmission parameter $\eta$ is high but strictly less than 1.0, simultaneous multipixel estimation using symmetric superpositions (like the HB states) is vastly superior to sequential single-parameter estimations. We have proven that the required quantum correlations survive standard scattering decoherence, keeping the overall resolution securely below the classical Standard Quantum Limit (SQL).

---
**Status**: This summary is locked as the primary conclusion for the final microscopy report. The independent loss model ($p=1$) is established as the physical baseline.
