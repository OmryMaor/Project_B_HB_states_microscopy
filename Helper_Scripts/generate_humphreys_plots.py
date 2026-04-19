"""
Self-contained Humphreys Fig 3 replication: QCRB vs N and QCRB vs d
for eta = 1.0, 0.85, 0.9, 0.95 under independent photon-loss noise.

Uses Kraus-operator approach for the loss channel (memory efficient).
All physics functions are inline — no imports from legacy Project A code.

Fig 3A: d=2 fixed, n=1..6  → N=3,6,9,12,15,18  (all dims ≤ 6,859)
Fig 3B: n=1 fixed, d=1..4  → N=2,3,4,5          (all dims ≤ 3,125)
"""
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend: no GUI windows

import numpy as np
import matplotlib.pyplot as plt
from itertools import product
from scipy.special import comb as binomial
from qutip import *
import gc
import os
import sys

# ---- Output folder ----
out_folder = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    "plots_output", "Humphreys_Replications"
)
os.makedirs(out_folder, exist_ok=True)

# ===================================================================
# PHYSICS FUNCTIONS
# ===================================================================

def generate_combinations(K, N):
    """Signal-mode Fock combinations with sum <= N."""
    return [c for c in product(range(N + 1), repeat=K) if sum(c) <= N]


def initialize_HB_state(n, d):
    """
    Construct the Holland-Burnett state HB(n, d+1).
    
    Parameters
    ----------
    n : int   – photons per mode in the input (total N = n*(d+1))
    d : int   – number of phase-sensing (signal) modes
    
    Returns
    -------
    coeffs : ndarray of complex  – Fock-basis coefficients (signal modes)
    combs  : list of tuples      – corresponding Fock combinations
    N      : int                 – total photon number
    """
    N = n * (d + 1)
    D = N + 1
    modes = d + 1

    omega = np.exp(2j * np.pi / modes)

    # Creation operators for each physical mode
    a_ops = []
    for i in range(modes):
        ops = [qeye(D)] * modes
        ops[i] = create(D)
        a_ops.append(tensor(ops))

    # DFT-transformed (beam-splitter network) mode operators
    b_ops = []
    for k in range(modes):
        b_k = 0
        for j in range(modes):
            b_k += (1.0 / np.sqrt(modes)) * (omega ** (j * k)) * a_ops[j]
        b_ops.append(b_k)

    # Build the state: put n photons in each b-mode
    state = tensor([basis(D, 0)] * modes)
    for k in range(modes):
        for _ in range(n):
            state = b_ops[k] * state
    state = state.unit()

    # Extract coefficients in the Fock basis of signal modes
    combinations = generate_combinations(d, N)
    coeffs = []
    for comb in combinations:
        n_0 = N - sum(comb)
        if n_0 < 0 or n_0 >= D:
            coeffs.append(0.0)
            continue
        target_ket = tensor([basis(D, n_0)] + [basis(D, ni) for ni in comb])
        coeffs.append(target_ket.overlap(state))

    return np.array(coeffs), combinations, N


def single_mode_loss_kraus(eta, D):
    """
    Kraus operators for a single-mode photon-loss channel with
    transmission probability eta, truncated to D Fock levels.
    
    E_k = sum_n sqrt(C(n,k)) * eta^((n-k)/2) * (1-eta)^(k/2) |n-k><n|
    for k = 0, 1, ..., D-1
    """
    kraus_ops = []
    for k in range(D):
        mat = np.zeros((D, D), dtype=complex)
        for n in range(k, D):
            mat[n - k, n] = np.sqrt(binomial(n, k, exact=True)) * \
                            eta ** ((n - k) / 2.0) * \
                            (1.0 - eta) ** (k / 2.0)
        kraus_ops.append(Qobj(mat))
    return kraus_ops


def apply_loss_kraus(rho_full, eta, d, D):
    """
    Apply independent photon-loss to each of the d signal modes (indices 1..d).
    Mode 0 (reference) is left lossless.
    
    Uses the Kraus representation — works entirely in the original
    (d+1)-mode Hilbert space without any ancilla.
    """
    modes = d + 1
    kraus_single = single_mode_loss_kraus(eta, D)

    rho = rho_full
    # Apply loss channel sequentially to each signal mode
    for sig_mode in range(1, modes):
        rho_new = rho * 0  # zero operator with same dims
        for E_k in kraus_single:
            # Embed E_k on signal mode sig_mode, identity everywhere else
            ops = [qeye(D)] * modes
            ops[sig_mode] = E_k
            E_full = tensor(ops)
            rho_new = rho_new + E_full * rho * E_full.dag()
        rho = rho_new
    return rho


def calculate_QFIM_mixed(rho, d, D):
    """
    Compute the d×d QFIM for the d signal-mode phases using the spectral
    decomposition of the (d+1)-mode density matrix rho.

    The generators are the number operators n_1, ..., n_d on signal modes.
    """
    modes = d + 1
    K = d

    # Number operators for signal modes 1..d
    generators = []
    for i in range(K):
        ops = [qeye(D)] * modes
        ops[i + 1] = num(D)      # mode 0 is ref, modes 1..d are signal
        generators.append(tensor(ops))

    vals, vecs = rho.eigenstates()

    # Change-of-basis matrix
    U = Qobj(np.hstack([v.full() for v in vecs]), dims=[rho.dims[0], rho.dims[0]])

    # Pre-compute derivative matrices in the eigenbasis
    diffs = vals[np.newaxis, :] - vals[:, np.newaxis]
    der_matrices = []
    for a in range(K):
        n_a_eigen = (U.dag() * generators[a] * U).full()
        der_matrices.append(-1j * diffs * n_a_eigen)

    dim = len(vals)
    qfim = np.zeros((K, K))
    for a in range(K):
        for b in range(a, K):
            s = 0.0
            for n in range(dim):
                for m in range(dim):
                    v_sum = vals[n] + vals[m]
                    if v_sum > 1e-14:
                        s += (2.0 / v_sum) * np.real(
                            der_matrices[a][n, m] * der_matrices[b][m, n]
                        )
            qfim[a, b] = s
            qfim[b, a] = s

    return qfim


def calculate_QFIM_pure(coeffs, combs, d):
    """Pure-state QFIM from Fock-basis coefficients (fast, no QuTiP)."""
    K = d
    n_exp = np.zeros(K)
    n_corr = np.zeros((K, K))
    for idx, c in enumerate(coeffs):
        p = np.abs(c) ** 2
        for a in range(K):
            n_exp[a] += p * combs[idx][a]
            for b in range(K):
                n_corr[a, b] += p * combs[idx][a] * combs[idx][b]
    qfim = np.zeros((K, K))
    for a in range(K):
        for b in range(K):
            qfim[a, b] = 4 * (n_corr[a, b] - n_exp[a] * n_exp[b])
    return qfim


# ===================================================================
# PARAMETERS
# ===================================================================
etas_to_test = [0.95, 0.90, 0.85]
colors = ['purple', 'orange', 'brown']

# Fig 3A: CRB vs N (d=2 fixed, n=1..6)
# d=2 keeps Hilbert dims small (max 19^3=6859) so ALL N points are computable.
# N = n*(d+1) = 3n  →  N = 3, 6, 9, 12, 15, 18
d_fixed = 2
n_values_a = [1, 2, 3, 4, 5, 6]
N_values_a = [n * (d_fixed + 1) for n in n_values_a]

# Fig 3B: CRB vs d (n=1 fixed, d=1..4)
n_fixed = 1
d_values_b = [1, 2, 3, 4]
N_values_b = [n_fixed * (d + 1) for d in d_values_b]

# Memory-safety: skip if (d+1)-mode Hilbert-space dim exceeds this
# Max expected dim: 19^3 = 6859 (n=6, d=2) — well within this cap
DIM_CAP = 8000

# ===================================================================
# 1.  PURE STATE BASELINES  (eta = 1.0)
# ===================================================================
print("=" * 60)
print("Computing pure-state (eta=1.0) baselines ...")
print("=" * 60)

crbs_pure_a = []
for n_val in n_values_a:
    c, cmbs, N_tot = initialize_HB_state(n_val, d_fixed)
    qf = calculate_QFIM_pure(c, cmbs, d_fixed)
    crbs_pure_a.append(np.real(np.trace(np.linalg.pinv(qf))))
    print(f"  HB({n_val},{d_fixed+1})  N={N_tot}  CRB={crbs_pure_a[-1]:.4f}")

crbs_pure_b = []
for d_val in d_values_b:
    c, cmbs, N_tot = initialize_HB_state(1, d_val)
    qf = calculate_QFIM_pure(c, cmbs, d_val)
    crbs_pure_b.append(np.real(np.trace(np.linalg.pinv(qf))))
    print(f"  HB(1,{d_val+1})  N={N_tot}  CRB={crbs_pure_b[-1]:.4f}")

# ===================================================================
# 2.  MIXED STATE SWEEP  (eta = 0.95, 0.9, 0.85)
# ===================================================================
print("\n" + "=" * 60)
print("Computing mixed-state QCRB (independent loss, Kraus) ...")
print("=" * 60)

# ---- Fig 3A data ----
mixed_a = {eta: [] for eta in etas_to_test}
valid_N = {eta: [] for eta in etas_to_test}

for eta in etas_to_test:
    for n_val in n_values_a:
        N_ph = n_val * (d_fixed + 1)
        D = N_ph + 1
        total_dim = D ** (d_fixed + 1)
        if total_dim > DIM_CAP:
            print(f"  [SKIP] eta={eta}, n={n_val}, N={N_ph}: dim={total_dim} > {DIM_CAP}")
            break  # higher n will be even larger
        print(f"  eta={eta}, n={n_val}, N={N_ph} (dim={total_dim}) ...", end=" ", flush=True)
        try:
            c, cmbs, _ = initialize_HB_state(n_val, d_fixed)
            # Reconstruct the full (d+1)-mode ket
            ket = 0
            for idx, c_val in enumerate(c):
                if np.abs(c_val) > 1e-12:
                    n_0 = N_ph - sum(cmbs[idx])
                    fock = tensor([basis(D, n_0)] + [basis(D, ni) for ni in cmbs[idx]])
                    ket = ket + c_val * fock
            rho_pure = ket2dm(ket.unit())

            rho_out = apply_loss_kraus(rho_pure, eta, d_fixed, D)
            qfim = calculate_QFIM_mixed(rho_out, d_fixed, D)
            crb = np.real(np.trace(np.linalg.pinv(qfim)))
            mixed_a[eta].append(crb)
            valid_N[eta].append(N_ph)
            print(f"CRB={crb:.4f}")
            gc.collect()
        except MemoryError:
            print("MemoryError — stopping this eta.")
            gc.collect()
            break
        except Exception as e:
            print(f"Error: {e}")
            import traceback; traceback.print_exc()
            break

# ---- Fig 3B data ----
mixed_b = {eta: [] for eta in etas_to_test}
valid_d = {eta: [] for eta in etas_to_test}

for eta in etas_to_test:
    for d_val in d_values_b:
        N_ph = n_fixed * (d_val + 1)
        D = N_ph + 1
        total_dim = D ** (d_val + 1)
        if total_dim > DIM_CAP:
            print(f"  [SKIP] eta={eta}, d={d_val}, N={N_ph}: dim={total_dim} > {DIM_CAP}")
            break
        print(f"  eta={eta}, d={d_val}, N={N_ph} (dim={total_dim}) ...", end=" ", flush=True)
        try:
            c, cmbs, _ = initialize_HB_state(n_fixed, d_val)
            ket = 0
            for idx, c_val in enumerate(c):
                if np.abs(c_val) > 1e-12:
                    n_0 = N_ph - sum(cmbs[idx])
                    fock = tensor([basis(D, n_0)] + [basis(D, ni) for ni in cmbs[idx]])
                    ket = ket + c_val * fock
            rho_pure = ket2dm(ket.unit())

            rho_out = apply_loss_kraus(rho_pure, eta, d_val, D)
            qfim = calculate_QFIM_mixed(rho_out, d_val, D)
            crb = np.real(np.trace(np.linalg.pinv(qfim)))
            mixed_b[eta].append(crb)
            valid_d[eta].append(d_val)
            print(f"CRB={crb:.4f}")
            gc.collect()
        except MemoryError:
            print("MemoryError — stopping this eta.")
            gc.collect()
            break
        except Exception as e:
            print(f"Error: {e}")
            import traceback; traceback.print_exc()
            break

# ===================================================================
# 3.  PLOTTING
# ===================================================================
print("\n" + "=" * 60)
print("Generating plots ...")
print("=" * 60)

# Ideal benchmarks
multipixel_a = [(d_fixed * (1 + np.sqrt(d_fixed))**2) / (4 * N**2) for N in N_values_a]
noon_a       = [(d_fixed**3) / N**2 for N in N_values_a]

multipixel_b = [(d * (1 + np.sqrt(d))**2) / (4 * N**2) for d, N in zip(d_values_b, N_values_b)]
noon_b       = [(d**3) / N**2 for d, N in zip(d_values_b, N_values_b)]

# ---- Fig 3A: QCRB vs N ----
fig, ax = plt.subplots(figsize=(8, 5))

# Pure state baseline
ax.plot(N_values_a, crbs_pure_a, 's-', color='green', markersize=8,
        label=r'HB(n,%d) $\eta=1.0$' % (d_fixed + 1))

# Mixed state curves
for i, eta in enumerate(etas_to_test):
    if valid_N[eta]:
        ax.plot(valid_N[eta], mixed_a[eta], 'o-', color=colors[i], markersize=8,
                label=r'HB(n,%d) $\eta=%s$' % (d_fixed + 1, eta))

# Ideal limits
ax.plot(N_values_a, multipixel_a, '--', color='red', lw=2,
        label=r'Optimal $|\psi_s\rangle$ (Ideal)')
ax.plot(N_values_a, noon_a, '--', color='blue', lw=2,
        label='N00N States (Ideal)')

ax.set_xlabel(f"Total photon number N = {d_fixed+1}n", fontsize=12)
ax.set_ylabel(r"Total variance $|\Delta\theta|^2$", fontsize=12)
ax.set_title(f"CRB vs N for d={d_fixed} under Noise", fontsize=14)
ax.set_xticks(N_values_a)
ax.grid(True, ls=':', alpha=0.6)
ax.legend(fontsize=10)

path_a = os.path.join(out_folder, "QCRB_vs_N_d2_noise.png")
fig.savefig(path_a, dpi=300, bbox_inches='tight')
plt.close(fig)
print(f"  Saved -> {path_a}")

# ---- Fig 3A (Zoom-in): QCRB vs N ----
fig, ax = plt.subplots(figsize=(8, 5))
N_zoom = N_values_a[2:]  # N=9, 12, 15, 18
ax.plot(N_zoom, crbs_pure_a[2:], 's-', color='green', markersize=8,
        label=r'HB(n,%d) $\eta=1.0$' % (d_fixed + 1))
for i, eta in enumerate(etas_to_test):
    if valid_N[eta]:
        ax.plot(N_zoom, mixed_a[eta][2:], 'o-', color=colors[i], markersize=8,
                label=r'HB(n,%d) $\eta=%s$' % (d_fixed + 1, eta))
ax.plot(N_zoom, multipixel_a[2:], '--', color='red', lw=2,
        label=r'Optimal $|\psi_s\rangle$ (Ideal)')
ax.plot(N_zoom, noon_a[2:], '--', color='blue', lw=2,
        label='N00N States (Ideal)')

ax.set_xlabel(f"Total photon number N = {d_fixed+1}n", fontsize=12)
ax.set_ylabel(r"Total variance $|\Delta\theta|^2$", fontsize=12)
ax.set_title(f"CRB vs N for d={d_fixed} (Zoom-in)", fontsize=14)
ax.set_xticks(N_zoom)
ax.set_ylim(0, 0.2)
ax.grid(True, ls=':', alpha=0.6)
ax.legend(fontsize=10)

path_zoom = os.path.join(out_folder, "QCRB_vs_N_d2_zoomin.png")
fig.savefig(path_zoom, dpi=300, bbox_inches='tight')
plt.close(fig)
print(f"  Saved -> {path_zoom}")

# ---- Fig 3B: QCRB vs d ----
fig, ax = plt.subplots(figsize=(8, 5))

ax.plot(d_values_b, crbs_pure_b, 's-', color='green', markersize=8,
        label=r'HB(1, d) $\eta=1.0$')

for i, eta in enumerate(etas_to_test):
    if valid_d[eta]:
        ax.plot(valid_d[eta], mixed_b[eta], 'o-', color=colors[i], markersize=8,
                label=r'HB(1, d) $\eta=%s$' % eta)

ax.plot(d_values_b, multipixel_b, '--', color='red', lw=2,
        label=r'Optimal $|\psi_s\rangle$ (Ideal)')
ax.plot(d_values_b, noon_b, '--', color='blue', lw=2,
        label='N00N States (Ideal)')

ax.set_xlabel("Number of phases (d)", fontsize=12)
ax.set_ylabel(r"Total variance $|\Delta\theta|^2$", fontsize=12)
ax.set_title(f"CRB vs d for n={n_fixed} under Noise", fontsize=14)
ax.set_xticks(d_values_b)
ax.grid(True, ls=':', alpha=0.6)
ax.legend(fontsize=10)

path_b = os.path.join(out_folder, "QCRB_vs_d_noise.png")
fig.savefig(path_b, dpi=300, bbox_inches='tight')
plt.close(fig)
print(f"  Saved -> {path_b}")

print("\nDone!")
