import nbformat
from nbformat.v4 import new_notebook, new_markdown_cell, new_code_cell

nb = new_notebook()

# Cell 0: Markdown
nb.cells.append(new_markdown_cell("""# HB States Metrology Optimization & Comparison
This notebook initializes Holland-Burnett (HB) states for N photons and d phases, computes their QFI and CRB, and natively compares them against previously obtained **Optimized States** under identical parameters (saved in `Project_A_files/optimized_states_common_eta.csv`)."""))

# Cell 1: Imports
imports = """import numpy as np
from itertools import product
from qutip import *
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.gridspec import GridSpec
import pandas as pd
import os
import ast
"""
nb.cells.append(new_code_cell(imports))

# Cell 2: Core Physics Helper
physics_helpers = r'''def generate_combinations(K, N):
    all_combinations = product(range(N + 1), repeat=K)
    return [comb for comb in all_combinations if sum(comb) <= N]

def diagonalize(rho):
    return rho.eigenstates()

def compute_rho_derivative(rho, vals, vecs, generators, a):
    n_a = generators[a]
    U = Qobj(np.hstack([v.full() for v in vecs]), dims=[rho.dims[0], rho.dims[0]])
    n_a_eigen = (U.dag() * n_a * U).full()
    diffs = vals[np.newaxis, :] - vals[:, np.newaxis]
    return -1j * diffs * n_a_eigen

def calculate_QFIM(rho, vals, vecs, generators):
    dim = len(vals)
    K = len(generators)
    qfim = np.zeros((K, K))
    derivatives_eigen = [compute_rho_derivative(rho, vals, vecs, generators, a) for a in range(K)]

    for a in range(K):
        for b in range(a, K):
            term_sum = 0
            for n in range(dim):
                for m in range(dim):
                    v_sum = vals[n] + vals[m]
                    if v_sum > 1e-14:
                        elem_a = derivatives_eigen[a][n, m]
                        elem_b = derivatives_eigen[b][m, n]
                        term_sum += (2.0 / v_sum) * np.real(elem_a * elem_b)
            qfim[a, b] = term_sum
            qfim[b, a] = term_sum
    return qfim

import math
def N_free(initial_signal_state, eta_vec, N, K):
    D = N + 1
    rho = ket2dm(initial_signal_state) if initial_signal_state.isket else initial_signal_state
    for i in range(K):
        eta = eta_vec[i]
        kraus_ops = []
        for k_loss in range(D):
            E_k = np.zeros((D, D))
            for n in range(k_loss, D):
                E_k[n-k_loss, n] = np.sqrt(math.comb(n, k_loss) * (eta**(n-k_loss)) * ((1-eta)**k_loss))
            ops = [qeye(D)] * K
            ops[i] = Qobj(E_k)
            kraus_ops.append(tensor(ops))
            
        rho_new = 0 * rho
        for E in kraus_ops:
            rho_new += E * rho * E.dag()
        rho = rho_new
    return rho

def N_int(initial_signal_state, eta_vec, N, K):
    D = N + 1
    a_e = tensor([qeye(D)] * K + [destroy(D)])
    h_tot = tensor([qeye(D)] * (K + 1)) * 0.0
    for i in range(K):
        theta_i = np.arccos(np.sqrt(eta_vec[i]))
        ops_s = [qeye(D)] * (K + 1); ops_s[i] = destroy(D)
        h_tot += theta_i * (tensor(ops_s).dag() * a_e - tensor(ops_s) * a_e.dag())
    V_total = h_tot.expm()
    rho_joint = ket2dm(tensor(initial_signal_state, basis(D, 0)))
    return (V_total * rho_joint * V_total.dag()).ptrace(list(range(K)))

def N_tot(initial_ket, eta_vec, p, N, K):
    return p * N_free(initial_ket, eta_vec, N, K) + (1.0 - p) * N_int(initial_ket, eta_vec, N, K)
'''
nb.cells.append(new_code_cell(physics_helpers))

# Cell 3: Initializations
initialization = r'''def initialize_HB_state(N, d, D=None):
    D = N + 1 if D is None else D
    modes = d + 1
    n_list = [N // modes] * modes
    for i in range(N % modes):
        n_list[i] += 1
        
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
        for _ in range(n_list[k]): state = b_ops[k] * state

    state = state.unit() 
    combinations = generate_combinations(d, N)
    basis_kets_K = [tensor([basis(D, n_signal) for n_signal in comb]) for comb in combinations]

    coeffs = []
    for comb in combinations:
        n_0 = N - sum(comb)
        target_ket = tensor([basis(D, n_0)] + [basis(D, n_i) for n_i in comb])
        coeffs.append(target_ket.overlap(state))
        
    return np.array(coeffs), basis_kets_K, combinations

def construct_ket(coeffs, basis_kets):
    return sum(coeffs[i] * basis_kets[i] for i in range(len(coeffs))).unit()
'''
nb.cells.append(new_code_cell(initialization))

# Cell 4: Evaluation and Overlay Plotting Logic
evals = r'''

def evaluate_state_metrics(initial_ket, N, K, eta_vec, p):
    rho_out = N_tot(initial_ket, eta_vec, p, N, K)
    vals, vecs = diagonalize(rho_out)
    
    D = N + 1
    generators = []
    for i in range(K):
        ops = [qeye(D)] * K
        ops[i] = num(D)
        generators.append(tensor(ops))
        
    qfim = calculate_QFIM(rho_out, vals, vecs, generators)
    max_qfi = np.trace(qfim)
    total_var = np.trace(np.linalg.pinv(qfim)) if np.linalg.det(qfim) != 0 else np.inf
    return max_qfi, total_var

def evaluate_grid(N, d, state_coeffs, basis_kets, eta_values, p_values):
    results = []
    initial_ket = construct_ket(state_coeffs, basis_kets)
    for eta in eta_values:
        for p in p_values:
            qfi, crb = evaluate_state_metrics(initial_ket, N, d, [eta]*d, p)
            results.append({"eta": eta, "p": p, "max_qfi": qfi, "total_var": crb})
            print(f"HB Computed eta={eta}, p={p:.1f} | QFI={qfi:.3f}, CRB={crb:.3f}")
    return pd.DataFrame(results)

def format_optimized_database(filepath):
    if not os.path.exists(filepath):
        print(f"File {filepath} not found.")
        return None
    df = pd.read_csv(filepath)
    def extract_eta(eta_str):
        clean = str(eta_str).replace('np.float64(', '').replace(')', '').replace('[', '').replace(']', '').split(',')
        return float(clean[0])
    df['eta'] = df['eta_vec'].apply(extract_eta)
    if 'total_variance' in df.columns:
        df.rename(columns={'total_variance': 'total_var'}, inplace=True)
    return df

def format_ket_latex(coeffs, combs):
    terms = []
    for c, comb in zip(coeffs, combs):
        if abs(c) > 0.01:
            c_cpx = complex(c)
            if abs(c_cpx.imag) < 1e-4:
                c_str = f"{c_cpx.real:.2f}"
            else:
                c_str = f"{c_cpx:.2f}"
            
            terms.append(f"{c_str}|{','.join(map(str, comb))}\\rangle")
    return " + ".join(terms).replace("+ -", "- ")

def plot_hb_vs_optimized(df_hb, df_opt, hb_coeffs, combs, N_target, K_target, eta_target, output_dir, file_prefix="HB_vs_Optimized"):
    plt.rcParams.update({'figure.facecolor': 'white', 'axes.facecolor': 'white', 'savefig.facecolor': 'white'})
    
    subset_hb = df_hb[np.isclose(df_hb['eta'], eta_target, atol=1e-3)].sort_values(by='p')
    if df_opt is not None:
        subset_opt = df_opt[np.isclose(df_opt['eta'], eta_target, atol=1e-3)].sort_values(by='p')
        valid_ps = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]
        subset_opt = subset_opt[pd.to_numeric(subset_opt['p']).apply(lambda x: any(np.isclose(x, v, atol=1e-3) for v in valid_ps))]
    else:
        subset_opt = pd.DataFrame()
        
    if subset_hb.empty and subset_opt.empty: return
        
    d, N = K_target, N_target
    noon_val = (d**3) / (N**2)
    computed_multipixel_crb = (d * (1 + np.sqrt(d))**2) / (4 * (N**2))
    
    fig, ax_qfi = plt.subplots(figsize=(10.5, 6.5))
    
    ket_labels_handles = []

    if not subset_opt.empty:
        colors_opt = plt.cm.viridis(np.linspace(0, 1, len(subset_opt)))
        for i, (_, row) in enumerate(subset_opt.iterrows()):
            lbl = r"CRB of $|\psi_{opt}(p)\rangle$" if i == 0 else ""
            scatter_h = ax_qfi.scatter(row['p'], row['total_var'], color=colors_opt[i], marker='*', s=300, edgecolors='black', zorder=4, alpha=0.9, label=lbl)
            
            c_str = str(row['coeffs']).split(',')
            cb_str = str(row['combinations']).split(';')
            try:
                c_arr = [float(x) for x in c_str if x.strip()]
                cb_arr = [tuple(map(int, x.split(','))) for x in cb_str if x.strip()]
                latex_ket = f"$p={row['p']:.1f}: {format_ket_latex(c_arr, cb_arr)}$"
                ket_labels_handles.append((scatter_h, latex_ket))
            except:
                pass


    if not subset_hb.empty:
        for i, (_, row) in enumerate(subset_hb.iterrows()):
            lbl = r"CRB of $|\psi_{HB}(p)\rangle$" if i == 0 else ""
            hb_h = ax_qfi.scatter(row['p'], row['total_var'], s=120, color='dodgerblue', edgecolors='black', zorder=5, label=lbl)

    ax_qfi.axhline(y=noon_val, color='red', linestyle='--', linewidth=2, zorder=2, label=f"Heisenberg Limit ($\\approx {noon_val:.2f}$)")
    ax_qfi.axhline(y=computed_multipixel_crb, color='dodgerblue', linestyle='-.', linewidth=2, zorder=2, label=f"Multipixel CRB ($\\approx {computed_multipixel_crb:.2f}$)")

    all_vars = list(subset_hb['total_var'].values) + ([noon_val, computed_multipixel_crb])
    if not subset_opt.empty: all_vars.extend(subset_opt['total_var'].values)
    v_min, v_max = min(all_vars), max(all_vars)
    padding = (v_max - v_min) * 0.15 if v_max > v_min else 0.2
    
    ax_qfi.set_ylim(v_min - padding * 0.5, v_max + padding * 3)

    ax_qfi.set_xlabel("Mixture Parameter $p$ (0 = Full Interaction, 1 = No Interaction)", fontsize=11)
    ax_qfi.set_ylabel(r"Minimized $\mathrm{Tr}(\mathcal{I}_{\boldsymbol{\theta}}^{-1})$", fontsize=11)
    ax_qfi.grid(True, linestyle=':', alpha=0.4)
    ax_qfi.set_title(f"State Comparison (HB vs Optimized CRB) for $\eta={eta_target}$", fontsize=13, pad=10)

    handles, labels = ax_qfi.get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    leg1 = ax_qfi.legend(by_label.values(), by_label.keys(), loc='upper left', fontsize=10, frameon=True, facecolor='white')
    ax_qfi.add_artist(leg1)
    
    l_hand = [h[0] for h in ket_labels_handles] if ket_labels_handles else []
    l_lab = [h[1] for h in ket_labels_handles] if ket_labels_handles else []
    
    if not subset_hb.empty:
        hb_ket_str = f"$|HB\\rangle: {format_ket_latex(hb_coeffs, combs)}$" 
        l_hand.append(hb_h)
        l_lab.append(hb_ket_str)
        
    if l_hand:
        title_str = r"$\mathbf{Optimized\ States\ (\eta=" + str(eta_target) + ")}$" if ket_labels_handles else r"$\mathbf{HB\ State\ (\eta=" + str(eta_target) + ")}$"
        leg2 = plt.legend(l_hand, l_lab, loc='center left', bbox_to_anchor=(1.02, 0.5), title=title_str, fontsize=9, frameon=True)
    
    plt.subplots_adjust(left=0.08, bottom=0.12, right=0.6, top=0.92)
    
    import os
    os.makedirs(output_dir, exist_ok=True)
    filename = f"{output_dir}/{file_prefix}_eta_{eta_target}.png"
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    print(f"Saved aesthetically matched scatter plot to {filename}")
    plt.close()

'''
nb.cells.append(new_code_cell(evals))

# Cell 5: Execution block
exec_block = r"""# == 1. Variables Definition ==
N_photons = 3
d_phases = 2
out_folder = "plots_output/HB_vs_Optimized"

hb_coeffs, basis_kets, combs = initialize_HB_state(N_photons, d_phases)

# == 2. Evaluation on Augmented Grid ==
p_grid = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]
# Requested explicit extended eta range constraints
eta_grid = [0.8, 0.9, 0.95, 1.0]

print("Evaluating HB State Matrix Profile...")
df_hb = evaluate_grid(N_photons, d_phases, hb_coeffs, basis_kets, eta_grid, p_grid)

# == 3. Load Project A Database ==
db_path = "Project_A_files/database_projectA_CRB.csv"
df_opt = format_optimized_database(db_path)

if df_opt is None:
    print("Warning: Project A database could not be loaded safely.")

# == 4. Generate Reproductions ==
for eta in eta_grid:
    plot_hb_vs_optimized(df_hb, df_opt, hb_coeffs, combs, N_photons, d_phases, eta, out_folder)
"""
nb.cells.append(new_code_cell(exec_block))

# Cell 6: Higher Dimensional HB Scans
exec_high_dims = r"""# == 5. Higher Dimension Sweeps ==
configs = [(4, 3), (5, 3), (6, 3)]
out_folder_high = "plots_output/HB_higher_dims"

import pandas as pd
import numpy as np

for (N_target, d_target) in configs:
    print(f"\n--- Evaluating HB State for N={N_target}, d={d_target} ---")
    hb_c, b_kets, cmbs = initialize_HB_state(N_target, d_target)
    
    # We dynamically re-evaluate for the grid
    df_hb_h = evaluate_grid(N_target, d_target, hb_c, b_kets, eta_grid, p_grid)
    for eta in eta_grid:
        plot_hb_vs_optimized(df_hb_h, None, hb_c, cmbs, N_target, d_target, eta, out_folder_high, file_prefix=f"HB_alone_N{N_target}_d{d_target}")
"""
nb.cells.append(new_code_cell(exec_high_dims))

# Write fully decoupled clean notebook 
with open("HB_QFI_loss.ipynb", "w", encoding='utf-8') as f:
    nbformat.write(nb, f)
print("Notebook HB_QFI_loss.ipynb Successfully Overwritten and Refactored!")
