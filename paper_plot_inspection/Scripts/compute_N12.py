"""
Compute N=12 (d=3) QCRB for eta=0.95, 0.9, 0.85 using the LOW-RANK trick.

The key insight: after applying Kraus loss operators to a pure state |ψ⟩,
the output ρ_out = Σ_k E_k|ψ⟩⟨ψ|E_k† has rank at most D^d (=2197 for N=12,d=3).
Instead of eigendecomposing the full 28,561×28,561 matrix, we eigendecompose
the 2,197×2,197 Gram matrix, which is ~170× smaller.
"""
import matplotlib
matplotlib.use('Agg')

import numpy as np
import matplotlib.pyplot as plt
from itertools import product
from scipy.special import comb as binomial
import gc
import os
import time

out_folder = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    "plots_output", "Humphreys_Replications"
)
os.makedirs(out_folder, exist_ok=True)

# ===================================================================
# PHYSICS FUNCTIONS — pure numpy, no QuTiP (for speed & memory)
# ===================================================================

def generate_combinations(K, N):
    return [c for c in product(range(N + 1), repeat=K) if sum(c) <= N]


def initialize_HB_state_numpy(n, d):
    """Build HB(n, d+1) state using numpy arrays (not QuTiP).
    Returns the state vector in the full (d+1)-mode Fock basis."""
    N = n * (d + 1)
    D = N + 1
    modes = d + 1
    total_dim = D ** modes

    omega = np.exp(2j * np.pi / modes)

    # DFT creation operator: b†_k = (1/sqrt(modes)) Σ_j ω^{jk} a†_j
    # Apply n times for each k-mode to the vacuum.
    # Build state vector directly.
    state = np.zeros(total_dim, dtype=complex)
    state[0] = 1.0  # vacuum

    def index_to_fock(idx, modes, D):
        """Convert linear index to Fock tuple (n_0, n_1, ..., n_{modes-1})."""
        fock = []
        for m in range(modes - 1, -1, -1):
            fock.append(idx // (D ** m))
            idx %= (D ** m)
        return tuple(fock)

    def fock_to_index(fock, D):
        """Convert Fock tuple to linear index."""
        idx = 0
        for m, n_m in enumerate(fock):
            idx += n_m * (D ** (len(fock) - 1 - m))
        return idx

    # Apply b†_k = (1/sqrt(modes)) Σ_j ω^{jk} a†_j, n times for each k
    for k in range(modes):
        for _ in range(n):
            new_state = np.zeros_like(state)
            for idx in range(total_dim):
                if abs(state[idx]) < 1e-15:
                    continue
                fock = list(index_to_fock(idx, modes, D))
                for j in range(modes):
                    if fock[j] + 1 < D:
                        coeff = (1.0 / np.sqrt(modes)) * (omega ** (j * k))
                        new_fock = fock.copy()
                        new_fock[j] += 1
                        new_idx = fock_to_index(new_fock, D)
                        new_state[new_idx] += coeff * np.sqrt(new_fock[j]) * state[idx]
            state = new_state

    # Normalize
    state /= np.linalg.norm(state)
    return state, D, N


def single_mode_kraus_matrices(eta, D):
    """Return list of D Kraus matrices (each D×D numpy array) for photon loss."""
    kraus = []
    for k in range(D):
        mat = np.zeros((D, D), dtype=complex)
        for n_val in range(k, D):
            mat[n_val - k, n_val] = np.sqrt(binomial(n_val, k, exact=True)) * \
                                     eta ** ((n_val - k) / 2.0) * \
                                     (1.0 - eta) ** (k / 2.0)
        kraus.append(mat)
    return kraus


def apply_single_mode_kraus_to_vectors(phi_vectors, kraus_mats, mode_idx, modes, D):
    """
    Apply single-mode Kraus operators on mode `mode_idx` to a set of state
    vectors. Each input vector |φ_i⟩ produces D output vectors E_k|φ_i⟩.

    Parameters
    ----------
    phi_vectors : list of 1D arrays (each of length D^modes)
    kraus_mats  : list of D×D numpy Kraus matrices
    mode_idx    : which mode to apply loss to
    modes       : total number of modes
    D           : local Hilbert space dimension

    Returns
    -------
    new_vectors : list of 1D arrays (D times as many as input)
    """
    new_vectors = []
    total_dim = D ** modes

    for phi in phi_vectors:
        # Reshape to tensor [D, D, ..., D]
        shape = [D] * modes
        tensor = phi.reshape(shape)

        for E_k in kraus_mats:
            # Apply E_k on axis `mode_idx`
            new_tensor = np.tensordot(E_k, tensor, axes=([1], [mode_idx]))
            # tensordot puts the E_k output axis first; move it to mode_idx
            new_tensor = np.moveaxis(new_tensor, 0, mode_idx)
            new_vectors.append(new_tensor.reshape(total_dim))

    return new_vectors


def compute_qcrb_lowrank(psi, eta, d, D):
    """
    Compute QCRB = Tr(QFIM^{-1}) for the HB state |ψ⟩ after independent
    photon loss with transmission η on each of the d signal modes.

    Uses the low-rank trick:
    1. Generate all output vectors |φ_k⟩ = E_k|ψ⟩ (composite Kraus)
    2. Form the Gram matrix G_{ij} = ⟨φ_i|φ_j⟩ (small: D^d × D^d)
    3. Eigendecompose G to get eigenvalues λ_i and coefficients
    4. Compute QFIM using only the non-zero-eigenvalue subspace
    """
    modes = d + 1
    total_dim = D ** modes
    t0 = time.time()

    # Step 1: Build all output vectors |φ_k⟩ by applying Kraus ops mode-by-mode
    kraus = single_mode_kraus_matrices(eta, D)
    print(f"    Building output vectors ...", end=" ", flush=True)

    phi_vectors = [psi.copy()]  # start with the pure state
    for sig_mode in range(1, modes):
        phi_vectors = apply_single_mode_kraus_to_vectors(
            phi_vectors, kraus, sig_mode, modes, D
        )
        gc.collect()

    n_vecs = len(phi_vectors)  # Should be D^d = 13^3 = 2197
    print(f"{n_vecs} vectors ({time.time()-t0:.1f}s)", flush=True)

    # Step 2: Form the matrix M (total_dim × n_vecs) where column i is |φ_i⟩
    # Then ρ_out = M M†
    # Instead of forming ρ_out, eigendecompose M† M (n_vecs × n_vecs)
    print(f"    Building Gram matrix ({n_vecs}×{n_vecs}) ...", end=" ", flush=True)
    t1 = time.time()

    # Stack into matrix: M has shape (total_dim, n_vecs)
    M = np.column_stack(phi_vectors)
    del phi_vectors
    gc.collect()

    # Gram matrix G = M† M
    G = M.conj().T @ M
    print(f"done ({time.time()-t1:.1f}s)", flush=True)

    # Step 3: Eigendecompose G
    print(f"    Eigendecomposing Gram matrix ...", end=" ", flush=True)
    t2 = time.time()
    eigvals, eigvecs_G = np.linalg.eigh(G)
    print(f"done ({time.time()-t2:.1f}s)", flush=True)

    # Keep only non-negligible eigenvalues
    mask = eigvals > 1e-12
    eigvals_r = eigvals[mask]
    eigvecs_G_r = eigvecs_G[:, mask]
    rank = len(eigvals_r)
    print(f"    Effective rank: {rank} / {n_vecs}", flush=True)

    # Step 4: Compute eigenvectors of ρ_out in the full space
    # v_i = M u_i / sqrt(λ_i)  where u_i are eigenvectors of G
    print(f"    Computing QFIM ...", end=" ", flush=True)
    t3 = time.time()

    V = M @ eigvecs_G_r  # shape (total_dim, rank)
    # Normalize: V[:, i] /= sqrt(eigvals_r[i])
    V = V / np.sqrt(eigvals_r)[np.newaxis, :]
    del M, G, eigvecs_G
    gc.collect()

    # Number operators for signal modes (modes 1..d)
    # For the QFIM we need ⟨v_i|n_a|v_j⟩ for each generator n_a
    # n_a acts on mode (a+1): n_a|n_0,...,n_a,...⟩ = n_{a+1} |n_0,...,n_a,...⟩
    # In the tensor index layout, this is just multiplication by the
    # occupation number of mode (a+1).

    K = d
    shape = [D] * modes

    # Precompute n_a matrix elements in the eigenbasis
    # n_a_eigen[i,j] = ⟨v_i|n_a|v_j⟩ = Σ_idx n_a(idx) * V[idx,i]* V[idx,j]
    # where n_a(idx) is the occupation number of mode a+1 at Fock index idx

    # Build the occupation-number vector for each signal mode
    occ_numbers = []
    for a in range(K):
        sig_mode = a + 1  # physical mode index (0=ref, 1..d=signal)
        occ = np.zeros(total_dim)
        for idx in range(total_dim):
            # fock[i] = (idx // D^(modes-1-i)) % D
            occ[idx] = (idx // (D ** (modes - 1 - sig_mode))) % D
        occ_numbers.append(occ)

    # Compute QFIM using vectorized operations
    # n_a_eigen = V† diag(occ_a) V
    qfim = np.zeros((K, K))

    # Weight matrix for QFIM
    sums = eigvals_r[:, None] + eigvals_r[None, :]
    weight = np.where(sums > 1e-14, 2.0 / sums, 0.0)
    diffs = eigvals_r[None, :] - eigvals_r[:, None]

    der_matrices = []
    for a in range(K):
        # n_a in eigenbasis: (rank × rank)
        # n_a_eigen = V† diag(occ_a) V = (occ_a[:, None] * V).T.conj() @ V
        occ_V = occ_numbers[a][:, np.newaxis] * V  # (total_dim, rank)
        n_a_eigen = V.conj().T @ occ_V  # (rank, rank)
        del occ_V
        der_matrices.append(-1j * diffs * n_a_eigen)

    for a in range(K):
        for b in range(a, K):
            val = np.sum(weight * np.real(der_matrices[a] * der_matrices[b].T))
            qfim[a, b] = val
            qfim[b, a] = val

    crb = np.real(np.trace(np.linalg.pinv(qfim)))
    print(f"done ({time.time()-t3:.1f}s)", flush=True)
    print(f"    Total time: {time.time()-t0:.1f}s", flush=True)

    return crb, qfim


# ===================================================================
# MAIN COMPUTATION
# ===================================================================
d_fixed = 3
etas = [0.95, 0.90, 0.85]
colors_eta = ['purple', 'orange', 'brown']

# Previously computed data
N_values_pure = [4, 8, 12]
crbs_pure_a = [0.7500, 0.2500, 0.1250]

d_values_b = [1, 2, 3, 4, 5, 6]
N_values_b = [2, 3, 4, 5, 6, 7]
crbs_pure_b = [0.2500, 0.5000, 0.7500, 1.0000, 1.2500, 1.5000]

known_a = {
    0.95: {'N': [4, 8], 'crb': [0.7907, 0.2706]},
    0.90: {'N': [4, 8], 'crb': [0.8385, 0.2953]},
    0.85: {'N': [4, 8], 'crb': [0.8951, 0.3252]},
}

known_b = {
    0.95: {'d': [1, 2, 3, 4], 'crb': [0.2635, 0.5270, 0.7907, 1.0543]},
    0.90: {'d': [1, 2, 3, 4], 'crb': [0.2793, 0.5585, 0.8385, 1.1184]},
    0.85: {'d': [1, 2, 3, 4], 'crb': [0.2980, 0.5956, 0.8951, 1.1946]},
}

# ===================================================================
# Build HB(3,4) state — N=12, d=3
# ===================================================================
print("=" * 60)
print("N=12, d=3: LOW-RANK approach")
print(f"Full Hilbert dim = 13^4 = 28,561")
print(f"Max output rank  = 13^3 = 2,197  (what we actually eigendecompose)")
print("=" * 60)

N_ph = 12
D = N_ph + 1  # 13

print("\nInitializing HB(3, 4) state ...", flush=True)
t_start = time.time()
psi, D_check, N_check = initialize_HB_state_numpy(3, d_fixed)
print(f"Done in {time.time()-t_start:.1f}s  (norm={np.linalg.norm(psi):.6f})")

for eta in etas:
    print(f"\n{'='*40}")
    print(f"  eta = {eta}")
    print(f"{'='*40}")
    crb, qfim = compute_qcrb_lowrank(psi, eta, d_fixed, D)
    print(f"  >>> CRB = {crb:.6f}")
    print(f"  >>> QFIM =\n{qfim}")

    known_a[eta]['N'].append(N_ph)
    known_a[eta]['crb'].append(crb)
    gc.collect()

# ===================================================================
# REGENERATE PLOTS
# ===================================================================
print("\n" + "=" * 60)
print("Regenerating final plots ...")
print("=" * 60)

# ---- Fig 3A: CRB vs N ----
fig, ax = plt.subplots(figsize=(8, 5))
ax.plot(N_values_pure, crbs_pure_a, 's-', color='green', markersize=8,
        label=r'HB(n, 4) $\eta=1.0$')
for i, eta in enumerate(etas):
    ax.plot(known_a[eta]['N'], known_a[eta]['crb'], 'o-', color=colors_eta[i],
            markersize=8, label=r'HB(n, 4) $\eta=%s$' % eta)
ax.plot(N_values_pure,
        [(d_fixed*(1+np.sqrt(d_fixed))**2)/(4*N**2) for N in N_values_pure],
        '--', color='red', lw=2, label=r'Optimal $|\psi_s\rangle$ (Ideal)')
ax.plot(N_values_pure,
        [(d_fixed**3)/N**2 for N in N_values_pure],
        '--', color='blue', lw=2, label='N00N States (Ideal)')
ax.set_xlabel('Number of photons (4n)', fontsize=12)
ax.set_ylabel(r'Total variance $|\Delta\theta|^2$', fontsize=12)
ax.set_title('CRB vs N for d=3 under Noise', fontsize=14)
ax.set_xticks(N_values_pure)
ax.grid(True, ls=':', alpha=0.6)
ax.legend()
fig.savefig(f'{out_folder}/Replication_Fig3a_Mixed_Loss.png', dpi=300, bbox_inches='tight')
plt.close(fig)
print(f"  Saved Replication_Fig3a_Mixed_Loss.png")

# ---- Fig 3B: CRB vs d (unchanged — d=5,6 infeasible) ----
fig, ax = plt.subplots(figsize=(8, 5))
ax.plot(d_values_b, crbs_pure_b, 's-', color='green', markersize=8,
        label=r'HB(1, d) $\eta=1.0$')
for i, eta in enumerate(etas):
    ax.plot(known_b[eta]['d'], known_b[eta]['crb'], 'o-', color=colors_eta[i],
            markersize=8, label=r'HB(1, d) $\eta=%s$' % eta)
ax.plot(d_values_b,
        [(d*(1+np.sqrt(d))**2)/(4*N**2) for d,N in zip(d_values_b, N_values_b)],
        '--', color='red', lw=2, label=r'Optimal $|\psi_s\rangle$ (Ideal)')
ax.plot(d_values_b,
        [d**3/N**2 for d,N in zip(d_values_b, N_values_b)],
        '--', color='blue', lw=2, label='N00N States (Ideal)')
ax.set_xlabel('Number of phases (d)', fontsize=12)
ax.set_ylabel(r'Total variance $|\Delta\theta|^2$', fontsize=12)
ax.set_title('CRB vs d for n=1 under Noise', fontsize=14)
ax.set_xticks(d_values_b)
ax.grid(True, ls=':', alpha=0.6)
ax.legend()
fig.savefig(f'{out_folder}/Replication_Fig3b_Mixed_Loss.png', dpi=300, bbox_inches='tight')
plt.close(fig)
print(f"  Saved Replication_Fig3b_Mixed_Loss.png")

print("\nAll done!")
