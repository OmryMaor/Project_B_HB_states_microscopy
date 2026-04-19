import numpy as np
import matplotlib.pyplot as plt
from itertools import product
from qutip import *
from Helper_Scripts.Archive_Code_Generation.fisher_A_source import N_tot, calculate_QFI_matrix_mixed, initialize_HB_state
import gc
import os

out_folder_high = "plots_output/Humphreys_Replications"
os.makedirs(out_folder_high, exist_ok=True)

etas_to_test = [0.95, 0.90, 0.85]
p_val = 1.0

d_fixed_a = 3
n_values_a = [1, 2, 3, 4]
N_values_a = [n * (d_fixed_a + 1) for n in n_values_a]

mixed_results_a = {eta: [] for eta in etas_to_test}
valid_N_a = {eta: [] for eta in etas_to_test}

print("Starting Mixed Noise Sweep for Fig 3A (CRB vs n)")
for eta in etas_to_test:
    for n, N_test in zip(n_values_a, N_values_a):
        print(f"  eta={eta}, N={N_test} (n={n})")
        try:
            hb_c, b_kets, cmbs = initialize_HB_state(N_test, d_fixed_a)
            D_loc = N_test + 1
            rho_pure = 0
            for idx, c in enumerate(hb_c):
                if np.abs(c) > 1e-8:
                    n_0 = N_test - sum(cmbs[idx])
                    ket = tensor([basis(D_loc, n_0)] + [basis(D_loc, n_i) for n_i in cmbs[idx]])
                    rho_pure += c * ket
            rho_pure = rho_pure * rho_pure.dag()
            
            rho_out = N_tot(p_val, eta, D_loc, d_fixed_a, rho_pure)
            qfim = calculate_QFI_matrix_mixed(rho_out, d_fixed_a, D_loc, cmbs, N_test)
            crb = np.trace(np.linalg.pinv(qfim))
            mixed_results_a[eta].append(crb)
            valid_N_a[eta].append(N_test)
        except MemoryError:
            print(f"    MemoryError at N={N_test}, gracefully capping limits.")
            gc.collect()
            break
        except Exception as e:
            print(f"    Error {e} at N={N_test}, capping limits.")
            break

plt.figure(figsize=(8, 5))
colors = ['purple', 'orange', 'brown']
for i, eta in enumerate(etas_to_test):
    if len(valid_N_a[eta]) > 0:
        plt.plot(valid_N_a[eta], mixed_results_a[eta], marker='o', linestyle='-', color=colors[i], label=f'HB(n,{d_fixed_a+1}) $\\eta={eta}$')

multipixel_crbs = [(d_fixed_a * (1 + np.sqrt(d_fixed_a))**2) / (4 * (n**2)) for n in N_values_a]
heisenberg_limits = [(d_fixed_a**3) / (n**2) for n in N_values_a]

plt.plot(N_values_a, multipixel_crbs, linestyle='--', color='red', linewidth=2, label=r'Optimal $|\psi_s\rangle$')
plt.plot(N_values_a, heisenberg_limits, linestyle='--', color='blue', linewidth=2, label='N00N States')
plt.xlabel(f"Number of photons ({d_fixed_a+1}n)", fontsize=12)
plt.ylabel(r"Total variance $|\Delta\theta|^2$", fontsize=12)
plt.title(f"CRB vs N for d={d_fixed_a} under Noise", fontsize=14)
plt.xticks(N_values_a)
plt.grid(True, linestyle=':', alpha=0.6)
plt.legend()
plt.savefig(f"{out_folder_high}/Replication_Fig3a_Mixed_Loss.png", dpi=300, bbox_inches='tight')
plt.close()

n_fixed_b = 1
d_values_b = [1, 2, 3, 4, 5, 6]
N_values_b = [n_fixed_b * (d + 1) for d in d_values_b]

mixed_results_b = {eta: [] for eta in etas_to_test}
valid_d_b = {eta: [] for eta in etas_to_test}

print("Starting Mixed Noise Sweep for Fig 3B (CRB vs d)")
for eta in etas_to_test:
    for d_test, N_t in zip(d_values_b, N_values_b):
        print(f"  eta={eta}, d={d_test}, N={N_t}")
        try:
            hb_c, b_kets, cmbs = initialize_HB_state(N_t, d_test)
            D_loc = N_t + 1
            rho_pure = 0
            for idx, c in enumerate(hb_c):
                if np.abs(c) > 1e-8:
                    n_0 = N_t - sum(cmbs[idx])
                    ket = tensor([basis(D_loc, n_0)] + [basis(D_loc, n_i) for n_i in cmbs[idx]])
                    rho_pure += c * ket
            rho_pure = rho_pure * rho_pure.dag()
            
            rho_out = N_tot(p_val, eta, D_loc, d_test, rho_pure)
            qfim = calculate_QFI_matrix_mixed(rho_out, d_test, D_loc, cmbs, N_t)
            crb = np.trace(np.linalg.pinv(qfim))
            mixed_results_b[eta].append(crb)
            valid_d_b[eta].append(d_test)
        except MemoryError:
            print(f"    MemoryError at d={d_test}, gracefully capping limits.")
            gc.collect()
            break
        except Exception as e:
            print(f"    Error {e} at d={d_test}, capping limits.")
            break

plt.figure(figsize=(8, 5))
for i, eta in enumerate(etas_to_test):
    if len(valid_d_b[eta]) > 0:
        plt.plot(valid_d_b[eta], mixed_results_b[eta], marker='o', linestyle='-', color=colors[i], label=f'HB(1, d) $\\eta={eta}$')

multipixel_d = [(d * (1 + np.sqrt(d))**2) / (4 * (N**2)) for d, N in zip(d_values_b, N_values_b)]
heisenberg_d = [(d**3) / (N**2) for d, N in zip(d_values_b, N_values_b)]

plt.plot(d_values_b, multipixel_d, linestyle='--', color='red', linewidth=2, label=r'Optimal $|\psi_s\rangle$')
plt.plot(d_values_b, heisenberg_d, linestyle='--', color='blue', linewidth=2, label='N00N States')
plt.xlabel("Number of phases (d)", fontsize=12)
plt.ylabel(r"Total variance $|\Delta\theta|^2$", fontsize=12)
plt.title(f"CRB vs d for n={n_fixed_b} under Noise", fontsize=14)
plt.xticks(d_values_b)
plt.grid(True, linestyle=':', alpha=0.6)
plt.legend()
plt.savefig(f"{out_folder_high}/Replication_Fig3b_Mixed_Loss.png", dpi=300, bbox_inches='tight')
plt.close()

print("Completed isolated Mixed state execution.")
