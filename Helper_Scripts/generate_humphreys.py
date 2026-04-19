import numpy as np
import matplotlib.pyplot as plt
from itertools import product
from qutip import *
import gc
import sys
import os

# Ensure successful import
sys.path.append(os.path.join(os.getcwd(), 'Helper_Scripts', 'Archive_Code_Generation'))
from temp_math import N_tot, calculate_QFI_matrix_mixed, initialize_HB_state
from test_pure import calculate_QFIM_pure

out_folder = "plots_output/Humphreys_Replications"
os.makedirs(out_folder, exist_ok=True)

# ---------------------------------------------------------
# FAST PURE STATE REPLICATION ETA=1.0
# ---------------------------------------------------------

d_fixed_a = 3
n_values_a = [1, 2, 3]
N_values_a = [n * (d_fixed_a + 1) for n in n_values_a]

crbs_pure_a = []
for N_test in N_values_a:
    c, cmbs = initialize_HB_state(N_test, d_fixed_a)
    qf = calculate_QFIM_pure(c, cmbs, d_fixed_a)
    crbs_pure_a.append(np.trace(np.linalg.pinv(qf)))

plt.figure(figsize=(8, 5))
plt.plot(N_values_a, crbs_pure_a, marker='o', linestyle='None', color='green', markersize=8, zorder=5, label=f'HB(n,{d_fixed_a+1})')
plt.plot(N_values_a, [(d_fixed_a * (1 + np.sqrt(d_fixed_a))**2) / (4 * (n**2)) for n in N_values_a], '--', color='red', lw=2, label=r'Optimal $|\psi_s\rangle$')
plt.plot(N_values_a, [(d_fixed_a**3) / (n**2) for n in N_values_a], '--', color='blue', lw=2, label='N00N States')
plt.xlabel(f"Number of photons ({d_fixed_a+1}n)", fontsize=12)
plt.ylabel(r"Total variance $|\Delta\theta|^2$", fontsize=12)
plt.title(f"QCRB vs N for d={d_fixed_a} ($\eta=1.0$)", fontsize=14)
plt.xticks(N_values_a)
plt.grid(True, linestyle=':', alpha=0.6)
plt.legend()
plt.savefig(f"{out_folder}/Replication_Fig3a_CRB_vs_N_pure.png", dpi=300, bbox_inches='tight')
plt.close()

n_fixed_b = 1
d_values_b = [1, 2, 3, 4]
N_values_b = [n_fixed_b * (d + 1) for d in d_values_b]

crbs_pure_b = []
for d_test, N_t in zip(d_values_b, N_values_b):
    c, cmbs = initialize_HB_state(N_t, d_test)
    qf = calculate_QFIM_pure(c, cmbs, d_test)
    crbs_pure_b.append(np.trace(np.linalg.pinv(qf)))

plt.figure(figsize=(8, 5))
plt.plot(d_values_b, crbs_pure_b, marker='o', linestyle='None', color='green', markersize=8, zorder=5, label=f'HB(1, d)')
plt.plot(d_values_b, [(d * (1 + np.sqrt(d))**2) / (4 * (N**2)) for d, N in zip(d_values_b, N_values_b)], '--', color='red', lw=2, label=r'Optimal $|\psi_s\rangle$')
plt.plot(d_values_b, [(d**3) / (N**2) for d, N in zip(d_values_b, N_values_b)], '--', color='blue', lw=2, label='N00N States')
plt.xlabel("Number of phases (d)", fontsize=12)
plt.ylabel(r"Total variance $|\Delta\theta|^2$", fontsize=12)
plt.title(f"QCRB vs d for n={n_fixed_b} ($\eta=1.0$)", fontsize=14)
plt.xticks(d_values_b)
plt.grid(True, linestyle=':', alpha=0.6)
plt.legend()
plt.savefig(f"{out_folder}/Replication_Fig3b_CRB_vs_d_pure.png", dpi=300, bbox_inches='tight')
plt.close()


# ---------------------------------------------------------
# MIXED STATE REPLICATION ETA=0.95, 0.9, 0.85
# ---------------------------------------------------------

etas_to_test = [0.95, 0.90, 0.85]
p_val = 1.0
colors = ['purple', 'orange', 'brown']

# FIG 3A Mixed
mixed_results_a = {eta: [] for eta in etas_to_test}
valid_N_a = {eta: [] for eta in etas_to_test}

for eta in etas_to_test:
    for n, N_test in zip(n_values_a, N_values_a):
        print(f"Mixed A: eta={eta}, N={N_test} (n={n})", flush=True)
        try:
            D_loc = N_test + 1
            if D_loc ** (d_fixed_a + 1) > 28600:
                print(f"  -> Dimension {D_loc**(d_fixed_a+1)} too large, graceful cap.", flush=True)
                break
            c, cmbs = initialize_HB_state(N_test, d_fixed_a)
            rho_pure = 0
            for idx, c_val in enumerate(c):
                if np.abs(c_val) > 1e-8:
                    n_0 = N_test - sum(cmbs[idx])
                    ket = tensor([basis(D_loc, n_0)] + [basis(D_loc, n_i) for n_i in cmbs[idx]])
                    rho_pure += c_val * ket
            if rho_pure != 0:
                rho_pure = rho_pure * rho_pure.dag()
            else:
                raise ValueError("Null state")
                
            rho_out = N_tot(p_val, eta, D_loc, d_fixed_a, rho_pure)
            qfim = calculate_QFI_matrix_mixed(rho_out, d_fixed_a, D_loc, cmbs, N_test)
            crb = np.trace(np.linalg.pinv(qfim))
            mixed_results_a[eta].append(crb)
            valid_N_a[eta].append(N_test)
        except MemoryError:
            print(f"  -> MemoryError at N={N_test}, graceful cap.")
            gc.collect()
            break
        except Exception as e:
            print(f"  -> Error {e} at N={N_test}, graceful cap.")
            break

plt.figure(figsize=(8, 5))
for i, eta in enumerate(etas_to_test):
    if len(valid_N_a[eta]) > 0:
        plt.plot(valid_N_a[eta], mixed_results_a[eta], marker='o', linestyle='-', color=colors[i], markersize=8, label=f'HB(n,{d_fixed_a+1}) $\eta={eta}$')
plt.plot(N_values_a, [(d_fixed_a * (1 + np.sqrt(d_fixed_a))**2) / (4 * (n**2)) for n in N_values_a], '--', color='red', lw=2, label=r'Optimal $|\psi_s\rangle$')
plt.plot(N_values_a, [(d_fixed_a**3) / (n**2) for n in N_values_a], '--', color='blue', lw=2, label='N00N States')
plt.xlabel(f"Number of photons ({d_fixed_a+1}n)", fontsize=12)
plt.ylabel(r"Total variance $|\Delta\theta|^2$", fontsize=12)
plt.title(f"QCRB vs N for d={d_fixed_a} under mixed noise", fontsize=14)
plt.xticks(N_values_a)
plt.grid(True, linestyle=':', alpha=0.6)
plt.legend()
plt.savefig(f"{out_folder}/Replication_Fig3a_Mixed_Loss.png", dpi=300, bbox_inches='tight')
plt.close()


# FIG 3B Mixed
mixed_results_b = {eta: [] for eta in etas_to_test}
valid_d_b = {eta: [] for eta in etas_to_test}

for eta in etas_to_test:
    for d_test, N_t in zip(d_values_b, N_values_b):
        print(f"Mixed B: eta={eta}, d={d_test}, N={N_t}", flush=True)
        try:
            D_loc = N_t + 1
            if D_loc ** (d_test + 1) > 28600:
                print(f"  -> Dimension {D_loc**(d_test+1)} too large, graceful cap.", flush=True)
                break
            c, cmbs = initialize_HB_state(N_t, d_test)
            rho_pure = 0
            for idx, c_val in enumerate(c):
                if np.abs(c_val) > 1e-8:
                    n_0 = N_t - sum(cmbs[idx])
                    ket = tensor([basis(D_loc, n_0)] + [basis(D_loc, n_i) for n_i in cmbs[idx]])
                    rho_pure += c_val * ket
            if rho_pure != 0:
                rho_pure = rho_pure * rho_pure.dag()
            else:
                raise ValueError("Null state")
                
            rho_out = N_tot(p_val, eta, D_loc, d_test, rho_pure)
            qfim = calculate_QFI_matrix_mixed(rho_out, d_test, D_loc, cmbs, N_t)
            crb = np.trace(np.linalg.pinv(qfim))
            mixed_results_b[eta].append(crb)
            valid_d_b[eta].append(d_test)
        except MemoryError:
            print(f"  -> MemoryError at d={d_test}, graceful cap.")
            gc.collect()
            break
        except Exception as e:
            print(f"  -> Error {e} at d={d_test}, graceful cap.")
            break

plt.figure(figsize=(8, 5))
for i, eta in enumerate(etas_to_test):
    if len(valid_d_b[eta]) > 0:
        plt.plot(valid_d_b[eta], mixed_results_b[eta], marker='o', linestyle='-', color=colors[i], markersize=8, label=f'HB(1, d) $\eta={eta}$')
plt.plot(d_values_b, [(d * (1 + np.sqrt(d))**2) / (4 * (N**2)) for d, N in zip(d_values_b, N_values_b)], '--', color='red', lw=2, label=r'Optimal $|\psi_s\rangle$')
plt.plot(d_values_b, [(d**3) / (N**2) for d, N in zip(d_values_b, N_values_b)], '--', color='blue', lw=2, label='N00N States')
plt.xlabel("Number of phases (d)", fontsize=12)
plt.ylabel(r"Total variance $|\Delta\theta|^2$", fontsize=12)
plt.title(f"QCRB vs d for n={n_fixed_b} under mixed noise", fontsize=14)
plt.xticks(d_values_b)
plt.grid(True, linestyle=':', alpha=0.6)
plt.legend()
plt.savefig(f"{out_folder}/Replication_Fig3b_Mixed_Loss.png", dpi=300, bbox_inches='tight')
plt.close()

print("All done!")
