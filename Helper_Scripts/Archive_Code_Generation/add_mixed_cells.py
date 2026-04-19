import json

notebook_path = "c:/Users/omrym/Technion/Semester9/Project B/HB_states_ProjectB/HB_QFI_loss.ipynb"

with open(notebook_path, 'r', encoding='utf-8') as f:
    nb = json.load(f)

# Find the start of the previous Fig 3 replication cells to keep them 
# and append a new one for mixed states.

new_cells = []

# Cell 4: Mixed states for Fig 3 replications (ETA=0.95, 0.9, 0.85)
code4 = """
# == REPLICATING FIG 3(A) & 3(B) UNDER MIXED NOISE (ETA = 0.95, 0.9, 0.85) ==
import gc

etas_to_test = [0.95, 0.90, 0.85]
p_val = 1.0 # Assuming non-correlated localized noise dominant for this scope

# --- FIG 3A: CRB vs N for fixed d=3 ---
d_fixed_a = 3
n_values_a = [1, 2, 3, 4]
N_values_a = [n * (d_fixed_a + 1) for n in n_values_a] # [4, 8, 12, 16]

mixed_results_a = {eta: [] for eta in etas_to_test}
valid_N_a = {eta: [] for eta in etas_to_test}

print("Starting Mixed Noise Sweep for Fig 3A (CRB vs n)")
for eta in etas_to_test:
    for n, N_test in zip(n_values_a, N_values_a):
        print(f"  eta={eta}, N={N_test} (n={n})")
        try:
            hb_c, b_kets, cmbs = initialize_HB_state(N_test, d_fixed_a)
            D_loc = N_test + 1
            
            # Reconstruct pure density matrix to feed into N_tot 
            # Note: memory heavy for N=12, 16
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
    plt.plot(valid_N_a[eta], mixed_results_a[eta], marker='o', linestyle='None', color=colors[i], label=f'HB(n,{d_fixed_a+1}) $\\eta={eta}$')

multipixel_crbs = [(d_fixed_a * (1 + np.sqrt(d_fixed_a))**2) / (4 * (n**2)) for n in N_values_a]
heisenberg_limits = [(d_fixed_a**3) / (n**2) for n in N_values_a]

plt.plot(N_values_a, multipixel_crbs, linestyle='--', color='red', linewidth=2, label=r'Optimal $|\psi_s\rangle$ (Ideal)')
plt.plot(N_values_a, heisenberg_limits, linestyle='--', color='blue', linewidth=2, label='N00N States (Ideal)')

plt.xlabel(f"Number of photons ({d_fixed_a+1}n)", fontsize=12)
plt.ylabel(r"Total variance $|\Delta\theta|^2$", fontsize=12)
plt.title(f"CRB vs N for d={d_fixed_a} under Noise", fontsize=14)
plt.xticks(N_values_a)
plt.grid(True, linestyle=':', alpha=0.6)
plt.legend()
plt.savefig(f"{out_folder_high}/Replication_Fig3a_Mixed_Loss.png", dpi=300, bbox_inches='tight')
pass # plt.show disabled

# --- FIG 3B: CRB vs d for fixed n=1 ---
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
    plt.plot(valid_d_b[eta], mixed_results_b[eta], marker='o', linestyle='None', color=colors[i], label=f'HB(1, d) $\\eta={eta}$')

multipixel_d = [(d * (1 + np.sqrt(d))**2) / (4 * (N**2)) for d, N in zip(d_values_b, N_values_b)]
heisenberg_d = [(d**3) / (N**2) for d, N in zip(d_values_b, N_values_b)]

plt.plot(d_values_b, multipixel_d, linestyle='--', color='red', linewidth=2, label=r'Optimal $|\psi_s\rangle$ (Ideal)')
plt.plot(d_values_b, heisenberg_d, linestyle='--', color='blue', linewidth=2, label='N00N States (Ideal)')

plt.xlabel("Number of phases (d)", fontsize=12)
plt.ylabel(r"Total variance $|\Delta\theta|^2$", fontsize=12)
plt.title(f"CRB vs d for n={n_fixed_b} under Noise", fontsize=14)
plt.xticks(d_values_b)
plt.grid(True, linestyle=':', alpha=0.6)
plt.legend()
plt.savefig(f"{out_folder_high}/Replication_Fig3b_Mixed_Loss.png", dpi=300, bbox_inches='tight')
pass # plt.show disabled

print("Completed Mixed state execution block!")
"""
new_cells.append({
    "cell_type": "code",
    "execution_count": None,
    "metadata": {},
    "outputs": [],
    "source": [line + "\n" for line in code4.split('\n')[:-1]] + [code4.split('\n')[-1]]
})

nb['cells'].extend(new_cells)

with open(notebook_path, 'w', encoding='utf-8') as f:
    json.dump(nb, f, indent=1)

print("Mixed simulation cells appended to notebook.")
