import numpy as np
import matplotlib.pyplot as plt
from qutip import *
from itertools import product
from scipy.special import comb as binomial
import gc
import os

# --- Setup ---
out_folder = "plots_output/Humphreys_Replications"
os.makedirs(out_folder, exist_ok=True)

# ===================================================================
# PHYSICS FUNCTIONS (Kraus-based)
# ===================================================================

def single_mode_loss_kraus(eta, D):
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
    modes = d + 1
    kraus_single = single_mode_loss_kraus(eta, D)
    rho = rho_full
    for sig_mode in range(1, modes):
        rho_new = rho * 0
        for E_k in kraus_single:
            ops = [qeye(D)] * modes
            ops[sig_mode] = E_k
            E_full = tensor(ops)
            rho_new = rho_new + E_full * rho * E_full.dag()
        rho = rho_new
    return rho

def calculate_QFIM_mixed(rho, d, D):
    modes = d + 1
    K = d
    generators = []
    for i in range(K):
        ops = [qeye(D)] * modes
        ops[i + 1] = num(D)
        generators.append(tensor(ops))
    vals, vecs = rho.eigenstates()
    U = Qobj(np.hstack([v.full() for v in vecs]), dims=[rho.dims[0], rho.dims[0]])
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
                        s += (2.0 / v_sum) * np.real(der_matrices[a][n, m] * der_matrices[b][m, n])
            qfim[a, b] = s
            qfim[b, a] = s
    return qfim

def initialize_HB_state(n, d):
    N = n * (d + 1)
    D = N + 1
    modes = d + 1
    omega = np.exp(2j * np.pi / modes)
    a_ops = []
    for i in range(modes):
        ops = [qeye(D)] * modes
        ops[i] = create(D)
        a_ops.append(tensor(ops))
    b_ops = []
    for k in range(modes):
        b_k = 0
        for j in range(modes):
            b_k += (1.0 / np.sqrt(modes)) * (omega ** (j * k)) * a_ops[j]
        b_ops.append(b_k)
    state = tensor([basis(D, 0)] * modes)
    for k in range(modes):
        for _ in range(n):
            state = b_ops[k] * state
    return state.unit(), N

def noisy_noon(N, d, eta):
    n_mode = N / d
    return (d**3 * (1 + eta**n_mode)) / (2 * (N**2) * (eta**n_mode))

# ===================================================================
# COMPUTATION (Matching QCRB_vs_N_d2_noise parameters)
# ===================================================================
etas = [0.95, 0.9, 0.85]
colors_eta = ['purple', 'orange', 'brown']

# FIG 3A: d=2, n=1..6
d_fixed = 2
n_vals_a = [1, 2, 3, 4, 5, 6]
hb_results_a = {eta: [] for eta in [1.0] + etas}
N_vals_a = []

print("Computing Plot A (d=2) ...")
for n in n_vals_a:
    psi, N_ph = initialize_HB_state(n, d_fixed)
    D_loc = N_ph + 1
    N_vals_a.append(N_ph)
    
    # Pure
    qf_pure = calculate_QFIM_mixed(ket2dm(psi), d_fixed, D_loc)
    hb_results_a[1.0].append(np.trace(np.linalg.pinv(qf_pure)))
    
    # Noisy
    rho_pure = ket2dm(psi)
    for eta in etas:
        rho_out = apply_loss_kraus(rho_pure, eta, d_fixed, D_loc)
        qf = calculate_QFIM_mixed(rho_out, d_fixed, D_loc)
        hb_results_a[eta].append(np.trace(np.linalg.pinv(qf)))
        gc.collect()
    print(f"  N={N_ph} done.")

# FIG 3B: n=1, d=1..4
n_fixed = 1
d_vals_b = [1, 2, 3, 4]
hb_results_b = {eta: [] for eta in [1.0] + etas}

print("\nComputing Plot B (n=1) ...")
for d in d_vals_b:
    psi, N_ph = initialize_HB_state(n_fixed, d)
    D_loc = N_ph + 1
    
    # Pure
    qf_pure = calculate_QFIM_mixed(ket2dm(psi), d, D_loc)
    hb_results_b[1.0].append(np.trace(np.linalg.pinv(qf_pure)))
    
    # Noisy
    rho_pure = ket2dm(psi)
    for eta in etas:
        rho_out = apply_loss_kraus(rho_pure, eta, d, D_loc)
        qf = calculate_QFIM_mixed(rho_out, d, D_loc)
        hb_results_b[eta].append(np.trace(np.linalg.pinv(qf)))
        gc.collect()
    print(f"  d={d} done.")

# ===================================================================
# PLOTTING
# ===================================================================

# PLOT A
plt.figure(figsize=(9, 6))
plt.plot(N_vals_a, hb_results_a[1.0], 's-', color='green', label=r'HB $\eta=1.0$')
for i, eta in enumerate(etas):
    plt.plot(N_vals_a, hb_results_a[eta], 'o-', color=colors_eta[i], label=r'HB $\eta=%s$' % eta)
    n00n_vals = [noisy_noon(N, d_fixed, eta) for N in N_vals_a]
    plt.plot(N_vals_a, n00n_vals, '--', color=colors_eta[i], alpha=0.7, label=r'N00N $\eta=%s$' % eta)

plt.xlabel(f'Total Photon Number N ({d_fixed+1}n)')
plt.ylabel('Total Variance')
plt.title(f'HB vs Noisy N00N: QCRB vs N (d={d_fixed})')
plt.grid(True, ls=':')
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()
plt.savefig(f"{out_folder}/QCRB_vs_N_d2_comparison_final.png", dpi=300)
plt.close()

# PLOT B
plt.figure(figsize=(9, 6))
plt.plot(d_vals_b, hb_results_b[1.0], 's-', color='green', label=r'HB $\eta=1.0$')
for i, eta in enumerate(etas):
    plt.plot(d_vals_b, hb_results_b[eta], 'o-', color=colors_eta[i], label=r'HB $\eta=%s$' % eta)
    n00n_vals = [noisy_noon(d+1, d, eta) for d in d_vals_b]
    plt.plot(d_vals_b, n00n_vals, '--', color=colors_eta[i], alpha=0.7, label=r'N00N $\eta=%s$' % eta)

plt.xlabel('Number of phases d (N = d+1)')
plt.ylabel('Total Variance')
plt.title('HB vs Noisy N00N: QCRB vs d (n=1)')
plt.grid(True, ls=':')
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()
plt.savefig(f"{out_folder}/QCRB_vs_d_comparison_final.png", dpi=300)
plt.close()

print(f"\nFinal consistency plots saved in {out_folder}")
