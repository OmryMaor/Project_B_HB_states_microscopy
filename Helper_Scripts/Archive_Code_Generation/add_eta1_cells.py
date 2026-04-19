import json
import os

notebook_path = "c:/Users/omrym/Technion/Semester9/Project B/HB_states_ProjectB/HB_QFI_loss.ipynb"

with open(notebook_path, 'r', encoding='utf-8') as f:
    nb = json.load(f)

new_cells = []

# Cell 1: Pure state evaluation method
code1 = """
# == FAST PURE STATE QFIM (ETA=1.0) ==

def calculate_QFIM_pure(coeffs, combs, d):
    # Mathematically exact and lightning fast QFIM calculation for pure states
    K = d
    qfim = np.zeros((K, K))
    n_exp = np.zeros(K)
    for idx, c in enumerate(coeffs):
        p = np.abs(c)**2
        for a in range(K):
            n_exp[a] += p * combs[idx][a]
            
    n_corr = np.zeros((K, K))
    for idx, c in enumerate(coeffs):
        p = np.abs(c)**2
        for a in range(K):
            for b in range(K):
                n_corr[a, b] += p * combs[idx][a] * combs[idx][b]
                
    for a in range(K):
        for b in range(K):
            qfim[a, b] = 4 * (n_corr[a, b] - n_exp[a] * n_exp[b])
    return qfim
"""
new_cells.append({
    "cell_type": "code",
    "execution_count": None,
    "metadata": {},
    "outputs": [],
    "source": [line + "\n" for line in code1.split('\n')[:-1]] + [code1.split('\n')[-1]]
})

# Cell 2: CRB vs N for d=4, eta=1 (Humphreys pure comparison)
code2 = """
# == NEW PLOT: CRB vs N (number of photons) for fixed d=4, eta=1.0 ==
# Using Pure state fast calculation to overcome large D scaling memory limits

d_fixed = 4
N_values_pure = [4, 8, 12, 16]
crb_N_pure_results = []
qfi_N_pure_results = []

print(f"Evaluating Fast Pure N-dependence for d={d_fixed}, eta=1.0...")
for N_test in N_values_pure:
    hb_c, b_kets, cmbs = initialize_HB_state(N_test, d_fixed)
    
    qf_mat = calculate_QFIM_pure(hb_c, cmbs, d_fixed)
    qfi = np.trace(qf_mat)
    crb = np.trace(np.linalg.pinv(qf_mat)) if np.linalg.det(qf_mat) != 0 else np.inf
    
    crb_N_pure_results.append(crb)
    qfi_N_pure_results.append(qfi)
    print(f"N={N_test} | QFI={qfi:.3f}, CRB={crb:.3f}")

plt.figure(figsize=(8, 5))
plt.plot(N_values_pure, crb_N_pure_results, marker='o', linestyle='None', color='darkorange', label=f'HB CRB ($\\eta=1.0$, Fast Pure)')

multipixel_crbs_pure = [(d_fixed * (1 + np.sqrt(d_fixed))**2) / (4 * (n**2)) for n in N_values_pure]
heisenberg_limits_pure = [(d_fixed**3) / (n**2) for n in N_values_pure]

plt.plot(N_values_pure, multipixel_crbs_pure, linestyle='-.', color='dodgerblue', label='Multipixel CRB (Ideal)')
plt.plot(N_values_pure, heisenberg_limits_pure, linestyle='--', color='red', label='Heisenberg Limit')

plt.xlabel("Number of Photons (N)", fontsize=12)
plt.ylabel(r"Minimum Variance (CRB)", fontsize=12)
plt.title(f"CRB vs N for d={d_fixed} ($\\eta=1.0, p=1.0$)", fontsize=14)
plt.grid(True, linestyle=':', alpha=0.6)
plt.legend()

filename_N_pure = f"{out_folder_high}/CRB_vs_N_d{d_fixed}_eta1.0_p1.0_scatter_pure.png"
plt.savefig(filename_N_pure, dpi=300, bbox_inches='tight')
print(f"Saved fast pure N-dependence plot to {filename_N_pure}")
pass # plt.show disabled
"""
new_cells.append({
    "cell_type": "code",
    "execution_count": None,
    "metadata": {},
    "outputs": [],
    "source": [line + "\n" for line in code2.split('\n')[:-1]] + [code2.split('\n')[-1]]
})

# Cell 3: CRB vs d for N=1, eta=1
code3 = """
# == NEW PLOT: CRB vs d for fixed N=1, eta=1.0 ==

N_fixed = 1
d_values_pure = [1, 2, 3, 4, 5, 6]
crb_d_pure_results = []

print(f"Evaluating Fast Pure d-dependence for N={N_fixed}, eta=1.0...")
for d_test in d_values_pure:
    hb_c, b_kets, cmbs = initialize_HB_state(N_fixed, d_test)
    
    qf_mat = calculate_QFIM_pure(hb_c, cmbs, d_test)
    # Using generalized inverse since determinant might be 0 for some edge cases like d=1
    crb = np.trace(np.linalg.pinv(qf_mat))
    
    crb_d_pure_results.append(crb)
    print(f"d={d_test} | CRB={crb:.3f}")

plt.figure(figsize=(8, 5))
plt.plot(d_values_pure, crb_d_pure_results, marker='s', linestyle='None', color='purple', label=f'HB CRB ($\\eta=1.0$, Fast Pure)')

multipixel_crbs_d_pure = [(d * (1 + np.sqrt(d))**2) / (4 * (N_fixed**2)) for d in d_values_pure]
heisenberg_limits_d_pure = [(d**3) / (N_fixed**2) for d in d_values_pure]

plt.plot(d_values_pure, multipixel_crbs_d_pure, linestyle='-.', color='dodgerblue', label='Multipixel CRB (Ideal)')
plt.plot(d_values_pure, heisenberg_limits_d_pure, linestyle='--', color='red', label='Heisenberg Limit')

plt.xlabel("Number of Pixels/Phases (d)", fontsize=12)
plt.ylabel(r"Minimum Variance (CRB)", fontsize=12)
plt.title(f"CRB vs d for N={N_fixed} ($\\eta=1.0, p=1.0$)", fontsize=14)
plt.xticks(d_values_pure)
plt.grid(True, linestyle=':', alpha=0.6)
plt.legend()

filename_d_pure = f"{out_folder_high}/CRB_vs_d_N{N_fixed}_eta1.0_p1.0_scatter_pure.png"
plt.savefig(filename_d_pure, dpi=300, bbox_inches='tight')
print(f"Saved fast pure d-dependence plot to {filename_d_pure}")
pass # plt.show disabled
"""
new_cells.append({
    "cell_type": "code",
    "execution_count": None,
    "metadata": {},
    "outputs": [],
    "source": [line + "\n" for line in code3.split('\n')[:-1]] + [code3.split('\n')[-1]]
})

nb['cells'].extend(new_cells)

with open(notebook_path, 'w', encoding='utf-8') as f:
    json.dump(nb, f, indent=1)

print("Cells for pure evaluation appended to notebook!")
