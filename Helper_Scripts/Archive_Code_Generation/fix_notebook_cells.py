import json

notebook_path = "c:/Users/omrym/Technion/Semester9/Project B/HB_states_ProjectB/HB_QFI_loss.ipynb"

with open(notebook_path, 'r', encoding='utf-8') as f:
    nb = json.load(f)

# Find where the pure functions were added
idx_start = -1
for i, cell in enumerate(nb['cells']):
    if cell['cell_type'] == 'code' and len(cell['source']) > 0 and '# == FAST PURE STATE' in cell['source'][0]:
        idx_start = i
        break

if idx_start != -1:
    # Remove the previously appended pure cells
    nb['cells'] = nb['cells'][:idx_start]

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

# Cell 2: CRB vs N for d=3, eta=1 (Humphreys pure comparison 3a)
code2 = """
# == REPLICATING FIG 3(A): CRB vs n for fixed d=3 (Labeled internally by authors as HB(n,4) typo) ==
# Total Photons N = 4n

d_fixed = 3
n_values = [1, 2, 3, 4]
N_values_pure = [n * (d_fixed + 1) for n in n_values] # [4, 8, 12, 16]

crb_N_pure_results = []

print(f"Evaluating Fast Pure N-dependence for d={d_fixed}, eta=1.0...")
for N_test in N_values_pure:
    hb_c, b_kets, cmbs = initialize_HB_state(N_test, d_fixed)
    qf_mat = calculate_QFIM_pure(hb_c, cmbs, d_fixed)
    crb = np.trace(np.linalg.pinv(qf_mat))
    crb_N_pure_results.append(crb)
    print(f"N={N_test} (n={N_test//4}) | CRB={crb:.3f}")

plt.figure(figsize=(8, 5))
# HB (green, points no line)
plt.plot(N_values_pure, crb_N_pure_results, marker='o', linestyle='None', color='green', label=f'HB(n,{d_fixed+1})')

multipixel_crbs_pure = [(d_fixed * (1 + np.sqrt(d_fixed))**2) / (4 * (n**2)) for n in N_values_pure]
heisenberg_limits_pure = [(d_fixed**3) / (n**2) for n in N_values_pure]

# Optimal Multipixel (red dashed)
plt.plot(N_values_pure, multipixel_crbs_pure, linestyle='--', color='red', linewidth=2, label=r'Optimal $|\psi_s\rangle$')
# Heisenberg Limit (N00N) (blue dashed)
plt.plot(N_values_pure, heisenberg_limits_pure, linestyle='--', color='blue', linewidth=2, label='N00N States')

plt.xlabel(f"Number of photons ({d_fixed+1}n)", fontsize=12)
plt.ylabel(r"Total variance $|\Delta\theta|^2$", fontsize=12)
plt.title(f"CRB vs N for d={d_fixed} (Replicating Fig 3a)", fontsize=14)
plt.xticks(N_values_pure)
plt.grid(True, linestyle=':', alpha=0.6)
plt.legend()

filename_N_pure = f"{out_folder_high}/Replication_Fig3a_CRB_vs_N_pure.png"
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

# Cell 3: CRB vs d for HB(1,d), eta=1 (Humphreys pure comparison 3b)
code3 = """
# == REPLICATING FIG 3(B): CRB vs d for fixed n=1 (Total Photons N = d+1) ==

n_fixed = 1
d_values_pure = [1, 2, 3, 4, 5, 6]
N_values_for_d = [n_fixed * (d + 1) for d in d_values_pure] # [2, 3, 4, 5, 6, 7]

crb_d_pure_results = []

print(f"Evaluating Fast Pure d-dependence for n={n_fixed}, eta=1.0...")
for d_test, N_test in zip(d_values_pure, N_values_for_d):
    hb_c, b_kets, cmbs = initialize_HB_state(N_test, d_test)
    qf_mat = calculate_QFIM_pure(hb_c, cmbs, d_test)
    crb = np.trace(np.linalg.pinv(qf_mat))
    crb_d_pure_results.append(crb)
    print(f"d={d_test}, N={N_test} | CRB={crb:.3f}")

plt.figure(figsize=(8, 5))
# HB (green, points no line)
plt.plot(d_values_pure, crb_d_pure_results, marker='o', linestyle='None', color='green', label=f'HB(1, d)')

multipixel_crbs_d_pure = [(d * (1 + np.sqrt(d))**2) / (4 * (N**2)) for d, N in zip(d_values_pure, N_values_for_d)]
heisenberg_limits_d_pure = [(d**3) / (N**2) for d, N in zip(d_values_pure, N_values_for_d)]

# Optimal Multipixel (red dashed)
plt.plot(d_values_pure, multipixel_crbs_d_pure, linestyle='--', color='red', linewidth=2, label=r'Optimal $|\psi_s\rangle$')
# Heisenberg Limit (N00N) (blue dashed)
plt.plot(d_values_pure, heisenberg_limits_d_pure, linestyle='--', color='blue', linewidth=2, label='N00N States')

plt.xlabel("Number of phases (d)", fontsize=12)
plt.ylabel(r"Total variance $|\Delta\theta|^2$", fontsize=12)
plt.title(f"CRB vs d for n={n_fixed} (Replicating Fig 3b)", fontsize=14)
plt.xticks(d_values_pure)
plt.grid(True, linestyle=':', alpha=0.6)
plt.legend()

filename_d_pure = f"{out_folder_high}/Replication_Fig3b_CRB_vs_d_pure.png"
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

print("Cells modified for proper Humprheys replication!")
