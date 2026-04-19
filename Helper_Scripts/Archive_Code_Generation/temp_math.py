import numpy as np
from qutip import *
from itertools import product
from scipy.optimize import minimize
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import pandas as pd
import os
import csv
import re
from matplotlib.gridspec import GridSpec

# The explicit database names for the 3 different runs
db_varying_eta = "optimized_states_varying_eta.csv"
db_common_opt = "optimized_states_common_eta.csv"
db_common_uni = "uniform_state.csv"

def generate_combinations(K, N):
    """
    Generate configurations where the sum of particles in K signal modes
    is less than or equal to N. The reference mode is implicit.
    """
    all_combinations = product(range(N + 1), repeat=K)
    # Corrected logic: sum(comb) <= N
    valid_combinations = [comb for comb in all_combinations if sum(comb) <= N]
    return valid_combinations

def diagonalize(rho):
    """
    Diagonalize the dm rho
    :param rho:
    :return: eigenvalues, eigenvectors
    """
    eigenvalues, eigenvectors = rho.eigenstates()
    return eigenvalues, eigenvectors

def compute_rho_derivative(rho, vals, vecs, generators, a):
    """
    Computes the matrix elements of the derivative ∂_a(ρ) = -1j[n_a, ρ] in the
    eigenbasis of ρ.

    Returns the matrix elements ⟨λ_n | ∂_a(ρ) | λ_m⟩ as a complex numpy array.

    :param rho: Density matrix (QuTiP Qobj)
    :param vals: Array of eigenvalues
    :param vecs: List of eigenvectors (QuTiP kets)
    :param generators: List of number operators
    :param a: Index of the parameter (mode) for the derivative
    :return: (np.ndarray) Matrix elements of the derivative in the eigenbasis
    """
    n_a = generators[a]

    # 1. Construct U, the diagonalizing matrix
    U = Qobj(np.hstack([v.full() for v in vecs]), dims=[rho.dims[0], rho.dims[0]])

    # 2. Transform the generator to the eigenbasis
    n_a_eigen = (U.dag() * n_a * U).full()

    # 3. Create the eigenvalue difference matrix (lambda_m - lambda_n)
    # diffs[n,m] = row[m] - col[n] = lambda_m - lambda_n
    diffs = vals[np.newaxis, :] - vals[:, np.newaxis]

    # 4. Return the result in the eigenbasis (No back-transformation)
    # This matrix contains the elements <n | del_a * rho | m>
    return -1j * diffs * n_a_eigen

def calculate_QFIM(rho, vals, vecs, generators):
    """
    Implements the Quantum Fisher Information Matrix (QFIM) using the spectral
    representation of the Symmetric Logarithmic Derivative (SLD).

    This function corresponds to Formula (11) in Albarelli et al. (2020),
    computing F_ab based on the matrix elements of ρ's derivatives in its
    own eigenbasis.
    https://iopscience.iop.org/article/10.1088/1751-8121/ab5d4d

    :param rho: Density matrix (Qobj)
    :param vals: Pre-calculated eigenvalues of rho
    :param vecs: Pre-calculated eigenvectors of rho
    :param generators: List of number operators for each phase mode
    :return: (np.ndarray) Symmetric KxK QFIM matrix
    """
    dim = len(vals)
    K = len(generators)
    qfim = np.zeros((K, K))


    # Pre-compute all derivative matrices in the eigenbasis
    derivatives_eigen = []
    for a in range(K):
        der_eigen = compute_rho_derivative(rho, vals, vecs, generators, a)
        derivatives_eigen.append(der_eigen)

    # Calculate QFIM elements. QFIM is symmetric, so iteration is on upper "triangle"
    for a in range(K):
        for b in range(a, K):
            term_sum = 0
            # Formula (11) double sum over eigenvalues
            for n in range(dim):
                for m in range(dim):
                    v_sum = vals[n] + vals[m]
                    if v_sum > 1e-14:
                        # Extract <n| del_a_rho |m> and <m| del_b_rho |n>
                        elem_a = derivatives_eigen[a][n, m]
                        elem_b = derivatives_eigen[b][m, n]

                        # Weighting 2 / (lambda_n + lambda_m)
                        term_sum += (2.0 / v_sum) * np.real(elem_a * elem_b)

            qfim[a, b] = term_sum
            qfim[b, a] = term_sum

    return qfim

def N_free(initial_signal_state, eta_vec, N, K):
    """
    Implements the independent noise channel (N_free) via Stinespring dilation.
    Now accepts eta_vec to allow different transmission coefficients per mode.
    """
    D = N + 1

    # 1. Global Unitary Initialization
    V_total = tensor([qeye(D)] * (2 * K))

    # 2. Sequential interaction application
    for i in range(K):
        # Calculate theta for this specific mode's eta
        theta_i = np.arccos(np.sqrt(eta_vec[i]))

        ops_s = [qeye(D)] * (2 * K); ops_s[i] = destroy(D)
        ops_e = [qeye(D)] * (2 * K); ops_e[K + i] = destroy(D)

        a_s = tensor(ops_s)
        a_e = tensor(ops_e)

        h_int = a_s.dag() * a_e - a_s * a_e.dag()
        V_total = (theta_i * h_int).expm() * V_total

    # 3. Initialize joint state
    env_vacuum = tensor([basis(D, 0)] * K)
    rho_joint = ket2dm(tensor(initial_signal_state, env_vacuum))

    rho_final_joint = V_total * rho_joint * V_total.dag()

    return rho_final_joint.ptrace(list(range(K)))

def N_int(initial_signal_state, eta_vec, N, K):
    """
    Implements N_int (Shared Environment).
    Signal: Modes 0 to K-1. Environment: Mode K.
    Total Modes: K + 1.
    Uses a vector of etas
    """
    D = N + 1
    a_e = tensor([qeye(D)] * K + [destroy(D)])

    h_tot = tensor([qeye(D)] * (K + 1)) * 0.0

    for i in range(K):
        # Calculate theta specific to this mode's transmission
        theta_i = np.arccos(np.sqrt(eta_vec[i]))

        ops_s = [qeye(D)] * (K + 1); ops_s[i] = destroy(D)
        a_s = tensor(ops_s)

        h_int_i = a_s.dag() * a_e - a_s * a_e.dag()

        # Add the properly weighted generator to the total Hamiltonian
        h_tot += theta_i * h_int_i

    # Exponentiate the properly summed generators
    V_total = h_tot.expm()

    rho_joint = ket2dm(tensor(initial_signal_state, basis(D, 0)))
    rho_final = V_total * rho_joint * V_total.dag()
    return rho_final.ptrace(list(range(K)))

def N_tot(initial_ket, eta_vec, p, N, K):
    """
    Implements the total noisy channel:
    rho_tot = p * N_free(rho_i) + (1 - p) * N_int(rho_i)

    Now accepts eta_vec to allow different transmission coefficients per mode.
    """
    # 1. Calculate the independent noise component (p probability)
    rho_free = N_free(initial_ket, eta_vec, N, K)

    # 2. Calculate the interacting noise component (1-p probability)
    rho_int = N_int(initial_ket, eta_vec, N, K)

    # 3. Sum the components to get the total mixed state
    rho_tot = p * rho_free + (1.0 - p) * rho_int

    return rho_tot

def optimize_initial_state(N, K, eta_vec, p):
    """
    Optimizes the initial probe state to maximize QFI under noise.
    Now accepts eta_vec for mode-specific losses and returns the total variance.
    """
    D = N + 1
    combinations = generate_combinations(K, N)
    num_coeffs = len(combinations)

    generators = []
    for i in range(K):
        ops = [qeye(D)] * K
        ops[i] = num(D)
        generators.append(tensor(ops))

    basis_kets = []
    for comb in combinations:
        state_modes = [basis(D, n) for n in comb]
        basis_kets.append(tensor(state_modes))

    def objective(coeffs):
        norm = np.linalg.norm(coeffs)
        if norm < 1e-10: return 0.0
        initial_ket = sum((coeffs[i] / norm) * basis_kets[i] for i in range(num_coeffs))

        # Updated to use eta_vec
        rho_out = N_tot(initial_ket, eta_vec, p, N, K)
        vals, vecs = diagonalize(rho_out)
        qfim = calculate_QFIM(rho_out, vals, vecs, generators)
        return -1.0 * np.trace(qfim)

    start_coeffs = np.array([1.0 / np.sqrt(num_coeffs)] * num_coeffs)
    cons = ({'type': 'eq', 'fun': lambda x: np.sum(x**2) - 1})
    bnds = [(0, 1) for _ in range(num_coeffs)]
    res = minimize(objective, start_coeffs, method='SLSQP', bounds=bnds, constraints=cons)

    # --- NEW: Evaluate the final metrics for the optimized state ---
    opt_coeffs = res.x
    norm = np.linalg.norm(opt_coeffs)
    initial_ket = sum((opt_coeffs[i] / norm) * basis_kets[i] for i in range(num_coeffs))

    rho_out = N_tot(initial_ket, eta_vec, p, N, K)
    vals, vecs = diagonalize(rho_out)
    qfim = calculate_QFIM(rho_out, vals, vecs, generators)

    max_qfi = np.trace(qfim)
    try:
        # Cramér-Rao Bound (Total Variance) is the trace of the inverse QFIM
        total_var = np.trace(np.linalg.pinv(qfim))
    except:
        total_var = np.inf

    return res.x, combinations, max_qfi, total_var

def format_ket_latex(coeffs, combinations, threshold=1e-2):
    """
    Helper to format kets using LaTeX notation with 2-digit precision.
    """
    terms = []
    for c, comb in zip(coeffs, combinations):
        if np.abs(c) > threshold:
            basis_str = f"|{','.join(map(str, comb))}\\rangle"
            terms.append(f"{c:.2f}{basis_str}")
    final_str = " + ".join(terms) if terms else "0"
    return f"${final_str}$"

def evaluate_uniform_state(N, K, eta_vec, p):
    """Evaluates the non-optimized uniform superposition state."""
    D = N + 1
    combinations = generate_combinations(K, N)
    num_coeffs = len(combinations)

    generators = []
    for i in range(K):
        ops = [qeye(D)] * K
        ops[i] = num(D)
        generators.append(tensor(ops))

    basis_kets = []
    for comb in combinations:
        state_modes = [basis(D, n) for n in comb]
        basis_kets.append(tensor(state_modes))

    coeffs = np.array([1.0 / np.sqrt(num_coeffs)] * num_coeffs)
    norm = np.linalg.norm(coeffs)
    initial_ket = sum((coeffs[i] / norm) * basis_kets[i] for i in range(num_coeffs))

    rho_out = N_tot(initial_ket, eta_vec, p, N, K)
    vals, vecs = diagonalize(rho_out)
    qfim = calculate_QFIM(rho_out, vals, vecs, generators)

    qfi = np.trace(qfim)
    try:
        total_var = np.trace(np.linalg.pinv(qfim))
    except:
        total_var = np.inf

    return coeffs, combinations, qfi, total_var

def reorganize_database_to_excel(input_file="optimized_states_database.csv", output_file="tidy_optimized_states.xlsx"):
    """
    Transforms the raw CSV into a formatted Excel table.
    Fock combinations are columns, coefficients are rounded to 6 decimal places.
    Updated for vector eta, total variance, and uniform state tracking.
    """
    if not os.path.exists(input_file):
        print(f"Error: '{input_file}' not found.")
        return

    # 1. Load and Parse
    df = pd.read_csv(input_file)
    tidy_data = []
    unique_combinations = set()

    for _, row in df.iterrows():
        entry = {
            'N': int(row['N']),
            'K': int(row['K']),
            'eta_vec': str(row['eta_vec']),  # Kept as string since it's a list
            'p': round(float(row['p']), 4),
            'max_qfi': round(float(row['max_qfi']), 6),
            'total_variance': round(float(row['total_variance']), 6) if not pd.isna(row['total_variance']) else 'inf',
            'is_uniform': bool(row['is_uniform'])
        }

        coeffs = [float(x) for x in str(row['coeffs']).split(',')]
        combs = [f"|{c}>" for c in str(row['combinations']).split(';')]

        for c, val in zip(combs, coeffs):
            entry[c] = round(val, 6)
            unique_combinations.add(c)
        tidy_data.append(entry)

    # 2. Structure DataFrame
    tidy_df = pd.DataFrame(tidy_data)
    metadata_cols = ['N', 'K', 'eta_vec', 'p', 'max_qfi', 'total_variance', 'is_uniform']
    basis_cols = sorted(list(unique_combinations))
    tidy_df = tidy_df.reindex(columns=metadata_cols + basis_cols).fillna(0.0)

    # 3. Save as Formatted Excel
    with pd.ExcelWriter(output_file, engine='xlsxwriter') as writer:
        tidy_df.to_excel(writer, index=False, sheet_name='Optimized Coefficients')

        workbook  = writer.book
        worksheet = writer.sheets['Optimized Coefficients']

        # Formatting: Freeze the first row and metadata columns (up to column 7)
        worksheet.freeze_panes(1, 7)

        # Apply a simple table style and adjust column widths
        header_format = workbook.add_format({'bold': True, 'bg_color': '#D7E4BC', 'border': 1})
        for col_num, value in enumerate(tidy_df.columns.values):
            worksheet.write(0, col_num, value, header_format)
            column_len = max(len(str(value)), 10)
            worksheet.set_column(col_num, col_num, column_len)

    print(f"--- Tidied Excel database saved as '{output_file}' ---")
    return tidy_df

# Execute
# tidy_db = reorganize_database_to_excel()

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os
import re
from matplotlib.gridspec import GridSpec

def format_ket_latex(coeffs, combinations, threshold=1e-2):
    terms = []
    for c, comb in zip(coeffs, combinations):
        if np.abs(c) > threshold:
            basis_str = f"|{','.join(map(str, comb))}\\rangle"
            terms.append(f"{c:.2f}{basis_str}")

    if not terms: return "$0$"
    chunks = [terms[i:i + 5] for i in range(0, len(terms), 5)]
    lines = []
    for i, chunk in enumerate(chunks):
        if i == 0: lines.append(f"${' + '.join(chunk)}$")
        else: lines.append(f"$+ {' + '.join(chunk)}$")
    return "\n".join(lines)

def # plot_from_database(N_target, K_target, eta_vec_target, filename, title_prefix="Optimized", save_path=None, forced_size=None, hide_legend=False, hide_benchmark_labels=False):
    if not os.path.exists(filename):
        return None

    plt.rcParams.update({
        'figure.facecolor': 'white', 'axes.facecolor': 'white', 'savefig.facecolor': 'white',
        'text.color': 'black', 'axes.labelcolor': 'black', 'xtick.color': 'black',
        'ytick.color': 'black', 'axes.edgecolor': 'black'
    })

    df = pd.read_csv(filename)
    target_scalar = float(eta_vec_target[0])

    def is_eta_match(s):
        clean_s = str(s).replace('np.float64(', '').replace(')', '')
        vals = re.findall(r"[-+]?(?:\d*\.\d+|\d+)(?:[eE][-+]?\d+)?", clean_s)
        if not vals: return False
        return np.isclose(float(vals[0]), target_scalar, atol=1e-3)

    p_targets = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]
    mask = (df['N'] == N_target) & (df['K'] == K_target) & \
           (df['eta_vec'].apply(is_eta_match)) & \
           (df['p'].apply(lambda x: any(np.isclose(x, pt, atol=0.01) for pt in p_targets)))

    subset = df[mask].sort_values(by='p')
    if subset.empty: return None

    current_size = forced_size if forced_size else 8.0
    fig = plt.figure(figsize=(current_size, current_size))
    fig.patch.set_facecolor('white')
    fig.patch.set_alpha(1.0)
    gs = GridSpec(2, 1, height_ratios=[4, 1], hspace=0.45)

    ax_qfi = fig.add_subplot(gs[0])
    ax_var = fig.add_subplot(gs[1])

    d, N = K_target, N_target
    sql_val = (d**2) / N
    noon_val = (d**3) / (N**2)
    optimal_hl_val = (d * (1 + np.sqrt(d))**2) / (4 * (N**2))

    # --- CALCULATE EXACT BOUNDS (LHS and RHS) ---
    computed_bound = sql_val # Safe fallback
    for _, row in subset.iterrows():
        if np.isclose(row['p'], 0.0, atol=0.01) and row['max_qfi'] > 0:
            computed_bound = (d**2) / row['max_qfi']
            break

    # Calculate Best Optimized State CRB
    best_qfi = subset['max_qfi'].max()
    best_crb = (d**2) / best_qfi

    lhs_bound = min(optimal_hl_val, computed_bound, best_crb)
    rhs_bound = max(sql_val, computed_bound, best_crb)

    margin = (rhs_bound - lhs_bound) * 0.1
    if margin == 0: margin = 0.1
    ax_var.set_xlim(lhs_bound - margin, rhs_bound + margin)
    ax_var.set_ylim(-1.0, 1.0)

    bar_height = 0.15
    ax_var.fill_between([lhs_bound, rhs_bound], -bar_height, bar_height, color='#E0E6ED', zorder=1)

    # --- MODIFIED: DECIMAL LABELS INSTEAD OF FRACTIONS ---
    benchmarks = [optimal_hl_val, noon_val, sql_val]
    benchmark_labels = [
        rf"$\mathbf{{Multipixel: \approx {optimal_hl_val:.2f}}}$",
        rf"$\mathbf{{Heisenberg: d^3/N^2 \approx {noon_val:.2f}}}$",
        rf"$\mathbf{{SQL: d^2/N \approx {sql_val:.2f}}}$"
    ]

    colors = plt.cm.viridis(np.linspace(0, 1, len(subset)))
    qfis = subset['max_qfi'].values

    # --- DRAW THE BENCHMARKS ---
    for i, (val, label) in enumerate(zip(benchmarks, benchmark_labels)):
        ax_var.plot([val, val], [-bar_height, bar_height], color='black', linewidth=2.5, zorder=2)

        if not hide_benchmark_labels:
            if i == 0: ax_var.text(val, bar_height + 0.1, label, color='black', fontsize=11, ha='left', va='bottom')
            elif i == 2: ax_var.text(val, bar_height + 0.1, label, color='black', fontsize=11, ha='right', va='bottom')
            else: ax_var.text(val, bar_height + 0.1, label, color='black', fontsize=11, ha='center', va='bottom')

    # --- PLOT THE DATA ---
    for i, (idx, row) in enumerate(subset.iterrows()):
        p_val = row['p']
        qfi_val = row['max_qfi']
        var_val = (d**2) / qfi_val

        coeffs = [float(x) for x in str(row['coeffs']).split(',')]
        combs = [tuple(map(int, c.split(','))) for c in str(row['combinations']).split(';')]

        latex_ket = format_ket_latex(coeffs, combs)
        legend_label = f"p={p_val:.2f}: {latex_ket}"

        ax_qfi.scatter(p_val, qfi_val, label=legend_label, s=120, color=colors[i], edgecolors='black', zorder=5)

        if np.isclose(p_val, 0.0, atol=0.01):
            ax_var.plot([var_val, var_val], [-bar_height, bar_height], color=colors[i], linewidth=8, solid_capstyle='butt', zorder=5)

    # --- MODIFIED: HIGHLIGHT BEST OPTIMIZED CRB ---
    ax_var.plot([best_crb, best_crb], [-bar_height, bar_height], color='darkorange', linewidth=4, linestyle='--', zorder=6)
    ax_var.text(best_crb, -bar_height - 0.1, rf"Best Opt: {best_crb:.2f}", color='darkorange', fontsize=11, ha='center', va='top', weight='bold')

    # --- MODIFIED: ADD HEISENBERG & MULTIPIXEL LIMIT QFI LINES ---
    hl_qfi_val = (d**2) / noon_val
    multipixel_qfi_val = (d**2) / optimal_hl_val

    ax_qfi.axhline(y=hl_qfi_val, color='red', linestyle='--', linewidth=2, zorder=2)
    ax_qfi.text(0.01, hl_qfi_val + 0.05, rf"Heisenberg Limit QFI ($\approx {hl_qfi_val:.2f}$)",
                color='red', fontsize=10, ha='left', va='bottom', weight='bold')

    ax_qfi.axhline(y=multipixel_qfi_val, color='dodgerblue', linestyle='-.', linewidth=2, zorder=2)
    ax_qfi.text(0.01, multipixel_qfi_val + 0.05, rf"Multipixel QFI ($\approx {multipixel_qfi_val:.2f}$)",
                color='dodgerblue', fontsize=10, ha='left', va='bottom', weight='bold')

    # Ensure all lines are inside the y-axis limits
    all_qfi_values = list(qfis) + [hl_qfi_val, multipixel_qfi_val]
    if len(all_qfi_values) > 0:
        q_min, q_max = np.min(all_qfi_values), np.max(all_qfi_values)
        if q_max == q_min: q_min -= 1; q_max += 1
        ax_qfi.set_ylim(q_min - 0.2*(q_max-q_min), q_max + 0.2*(q_max-q_min))

    ax_qfi.set_xlabel(r"Mixture Parameter $p$ (0 = Full Interaction, 1 = No Interaction)", fontsize=11)
    ax_qfi.set_ylabel(r"Maximized $\mathrm{Tr}(\mathcal{I}_{\boldsymbol{\theta}})$", fontsize=11)
    ax_qfi.grid(True, linestyle=':', alpha=0.4)

    ax_var.set_yticks([]); ax_var.spines['left'].set_visible(False)
    ax_var.spines['right'].set_visible(False); ax_var.spines['top'].set_visible(False)
    ax_var.spines['bottom'].set_visible(False); ax_var.set_xticks([])
    ax_var.set_title(r"Cramer-Rao Bound of State at $p=0$", fontsize=13, pad=5)

    if not hide_legend:
        handles, labels = ax_qfi.get_legend_handles_labels()
        plt.subplots_adjust(left=0.15, bottom=0.35, right=0.95, top=0.95)
        fig.legend(handles, labels, title=rf"{title_prefix} States $|\psi\rangle$",
                   loc='upper center', bbox_to_anchor=(0.5, 0.33),
                   fontsize=9, frameon=True, facecolor='white', ncol=2, title_fontsize='10')
    else:
        plt.subplots_adjust(left=0.15, bottom=0.15, right=0.95, top=0.95)

    if forced_size is None:
        fig.canvas.draw()
        bbox = fig.get_tightbbox(fig.canvas.get_renderer())
        plt.close(fig)
        return bbox.width, bbox.height
    else:
        if save_path:
            plt.savefig(save_path, dpi=300, facecolor='white', transparent=False)
            print(f"  -> Saved perfectly sized figure to {save_path}")
        pass # plt.show disabled
        return None

# ==========================================================
# 2-PASS PLOTTING EXECUTION
# ==========================================================
db_varying_eta = "optimized_states_varying_eta.csv"
db_common_opt = "optimized_states_common_eta.csv"
db_common_uni = "uniform_state.csv"
output_dir = "plots_output"
os.makedirs(output_dir, exist_ok=True)

print("--- PASS 1: Scanning for Maximum Figure Dimensions ---")
dim_opt = []
dim_uni = []

varying_etas_to_run = [[0.5, 0.5], [1.0, 0.5]]

for eta_vec_var in varying_etas_to_run:
    res = # plot_from_database(N_target=3, K_target=2, eta_vec_target=eta_vec_var, filename=db_varying_eta)
    if res: dim_opt.append(res)

for e in [0.1, 0.3, 0.6, 0.8]:
    res = # plot_from_database(N_target=3, K_target=2, eta_vec_target=[e, e], filename=db_common_opt)
    if res: dim_opt.append(res)

for e in [0.3, 0.8]:
    hide_bench = (e == 0.3)
    res = # plot_from_database(N_target=3, K_target=2, eta_vec_target=[e, e], filename=db_common_uni, hide_legend=True, hide_benchmark_labels=hide_bench)
    if res: dim_uni.append(res)

g_size_opt = max(max([d[0] for d in dim_opt]), max([d[1] for d in dim_opt])) + 0.5 if dim_opt else 9.0
g_size_uni = max(max([d[0] for d in dim_uni]), max([d[1] for d in dim_uni])) + 0.5 if dim_uni else 8.0

print("\n--- PASS 2: Generating Final Uncropped Plots ---")

for eta_vec_var in varying_etas_to_run:
    print(f"Plotting Varying Eta {eta_vec_var}...")
    # plot_from_database(N_target=3, K_target=2, eta_vec_target=eta_vec_var, filename=db_varying_eta,
                       title_prefix="Optimized Varying",
                       save_path=f"{output_dir}/varying_eta_{eta_vec_var[0]}_{eta_vec_var[1]}_custom_axi.png",
                       forced_size=g_size_opt)

for e in [0.3, 0.8]:
    print(f"Plotting Uniform State for Eta {e}...")
    hide_bench = (e == 0.3)
    # plot_from_database(N_target=3, K_target=2, eta_vec_target=[e, e], filename=db_common_uni,
                       title_prefix="Uniform",
                       save_path=f"{output_dir}/common_uni_eta_{e}_custom_axi.png",
                       forced_size=g_size_uni, hide_legend=True, hide_benchmark_labels=hide_bench)

for e in [0.1, 0.3, 0.6, 0.8]:
    print(f"Plotting Common Eta {e} (Optimized)...")
    # plot_from_database(N_target=3, K_target=2, eta_vec_target=[e, e], filename=db_common_opt,
                       title_prefix="Optimized Common",
                       save_path=f"{output_dir}/common_opt_eta_{e}_custom_axi.png",
                       forced_size=g_size_opt)

print("\nFinished generating split-logic 1D plots with min/max CRB bounds!")

import numpy as np
import csv
import os
import pandas as pd

N_param, K_param = 3, 2
db_varying_eta = "optimized_states_varying_eta.csv"
db_common_opt = "optimized_states_common_eta.csv"
db_common_uni = "uniform_state.csv"

def setup_csv_if_needed(filename):
    if not os.path.isfile(filename):
        with open(filename, mode='w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(['N', 'K', 'eta_vec', 'p', 'max_qfi', 'total_variance', 'coeffs', 'combinations'])

def get_completed_runs(filename):
    if not os.path.isfile(filename):
        return []
    df = pd.read_csv(filename)
    return list(zip(df['eta_vec'].astype(str), df['p']))

def is_computed(completed_list, eta_vec, p_target, tol=0.01):
    # FIX: Convert the target vector to pure Python floats before turning into a string!
    clean_target = str([float(e) for e in eta_vec])

    for ext_str, p_val in completed_list:
        clean_ext = ext_str.replace('np.float64(', '').replace(')', '')
        if clean_ext == clean_target and abs(float(p_val) - p_target) <= tol:
            return True
    return False

setup_csv_if_needed(db_varying_eta)
setup_csv_if_needed(db_common_opt)
setup_csv_if_needed(db_common_uni)

# ==========================================================
# RUN 1: Varying Eta Loop
# ==========================================================
varying_etas_to_run = [[0.5, 0.5], [1.0, 0.5]]
p_targets_var = [round(x, 1) for x in np.linspace(0.0, 1.0, 11)]
completed_var = get_completed_runs(db_varying_eta)

with open(db_varying_eta, mode='a', newline='') as f:
    writer = csv.writer(f)
    for eta_vec_var in varying_etas_to_run:
        print(f"\n--- Starting RUN 1: Varying Eta {eta_vec_var} ---")
        for p in p_targets_var:
            if is_computed(completed_var, eta_vec_var, p):
                print(f"  Skipping p = {p:.2f} (Already computed)")
                continue

            print(f"  Optimizing for p = {p:.2f}...")
            opt_coeffs, combs, max_qfi, opt_var = optimize_initial_state(N_param, K_param, eta_vec_var, p)
            writer.writerow([N_param, K_param, str(eta_vec_var), p, max_qfi, opt_var,
                             ",".join(map(str, opt_coeffs)),
                             ";".join([",".join(map(str, c)) for c in combs])])

# ==========================================================
# RUN 2 & 3: Common Eta Grid
# ==========================================================
print("\n--- Starting RUNS 2 & 3: Common Eta Grid (0.1 to 1.0) ---")
grid_values = [round(x, 1) for x in np.linspace(0.1, 1.0, 10)]
p_grid_values = [round(x, 1) for x in np.linspace(0.0, 1.0, 11)]

completed_opt = get_completed_runs(db_common_opt)
completed_uni = get_completed_runs(db_common_uni)

with open(db_common_opt, mode='a', newline='') as f_opt, \
     open(db_common_uni, mode='a', newline='') as f_uni:

    writer_opt = csv.writer(f_opt)
    writer_uni = csv.writer(f_uni)

    for eta_val in grid_values:
        # Convert np.float64 explicitly to standard float to avoid string issues
        current_eta_vec = [float(eta_val), float(eta_val)]
        print(f"\n> Sweeping common eta = {eta_val}")

        for p_val in p_grid_values:
            # Check and Run 2: Optimized State
            if is_computed(completed_opt, current_eta_vec, p_val):
                 print(f"    Skipping Optimized p = {p_val:.2f} (Already computed)")
            else:
                 print(f"    Optimizing p = {p_val:.2f}...")
                 opt_coeffs, combs, max_qfi, opt_var = optimize_initial_state(N_param, K_param, current_eta_vec, p_val)
                 writer_opt.writerow([N_param, K_param, str(current_eta_vec), p_val, max_qfi, opt_var,
                                      ",".join(map(str, opt_coeffs)),
                                      ";".join([",".join(map(str, c)) for c in combs])])

            # Check and Run 3: Uniform State
            if is_computed(completed_uni, current_eta_vec, p_val):
                 print(f"    Skipping Uniform p = {p_val:.2f} (Already computed)")
            else:
                 print(f"    Evaluating Uniform p = {p_val:.2f}...")
                 uni_coeffs, _, uni_qfi, uni_var = evaluate_uniform_state(N_param, K_param, current_eta_vec, p_val)
                 writer_uni.writerow([N_param, K_param, str(current_eta_vec), p_val, uni_qfi, uni_var,
                                      ",".join(map(str, uni_coeffs)),
                                      ";".join([",".join(map(str, c)) for c in combs])])

print("\n--- All experiments complete and appended to the databases! ---")

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import re
import os

def # plot_optimal_state_histogram(N_target=3, K_target=2, eta_vec_target=[0.6, 0.6], p_target=0.8, filename="optimized_states_common_eta.csv"):

    # FORCE PYCHARM MATPLOTLIB SETTINGS TO WHITE
    plt.rcParams.update({
        'figure.facecolor': 'white',
        'axes.facecolor': 'white',
        'savefig.facecolor': 'white',
        'text.color': 'black',
        'axes.labelcolor': 'black',
        'xtick.color': 'black',
        'ytick.color': 'black',
        'axes.edgecolor': 'black'
    })

    if not os.path.exists(filename):
        print(f"Database file '{filename}' not found.")
        return

    df = pd.read_csv(filename)

    # Robust matching (same as the 1D plots)
    target_scalar = float(eta_vec_target[0])
    def is_eta_match(s):
        clean_s = str(s).replace('np.float64(', '').replace(')', '')
        vals = re.findall(r"[-+]?(?:\d*\.\d+|\d+)(?:[eE][-+]?\d+)?", clean_s)
        if not vals: return False
        return np.isclose(float(vals[0]), target_scalar, atol=1e-3)

    mask = (df['N'] == N_target) & \
           (df['K'] == K_target) & \
           (df['eta_vec'].apply(is_eta_match)) & \
           (df['p'].apply(lambda x: np.isclose(float(x), p_target, atol=0.01)))

    subset = df[mask]

    if subset.empty:
        print(f"State for N={N_target}, K={K_target}, eta={eta_vec_target}, p={p_target} not found in database!")
        return

    # Extract coefficients and compute probability
    row = subset.iloc[0]
    coeffs = [float(x) for x in str(row['coeffs']).split(',')]
    combs = [tuple(map(int, c.split(','))) for c in str(row['combinations']).split(';')]

    probs = [c**2 for c in coeffs]
    labels = [f"$|{','.join(map(str, comb))}\\rangle$" for comb in combs]

    # Create the Plot (10x6 ratio to give the 10 labels room to breathe)
    fig = plt.figure(figsize=(10, 6))
    fig.patch.set_facecolor('white')
    fig.patch.set_alpha(1.0)

    ax = fig.add_subplot(111)

    # Magenta bars with a clean black edge
    bars = ax.bar(labels, probs, color='magenta', edgecolor='black', linewidth=1.5, zorder=3)

    # Formatting
    ax.grid(axis='y', linestyle='--', alpha=0.7, zorder=0)
    ax.set_xlabel("Basis States", fontsize=13, labelpad=10)
    ax.set_ylabel(r"Probability $|c_i|^2$", fontsize=13, labelpad=10)

    eta_str = eta_vec_target[0] if eta_vec_target[0] == eta_vec_target[1] else str(eta_vec_target)
    ax.set_title(f"Optimal State Probability Distribution\n(N={N_target}, K={K_target}, $\\eta$={eta_str}, p={p_target})", fontsize=15, pad=15)

    ax.set_ylim(0, max(probs) * 1.15) # Add a little headroom above the tallest bar

    # Add exact values on top of the bars
    for bar in bars:
        height = bar.get_height()
        if height > 0.005: # Only label visible bars (ignores exact 0s)
            ax.text(bar.get_x() + bar.get_width()/2., height + 0.01,
                    f'{height:.3f}', ha='center', va='bottom', fontsize=10)

    # Rotate labels so they don't overlap
    plt.xticks(rotation=45, ha='right', fontsize=11)

    # Fixed margins so the saved file size is consistent
    plt.subplots_adjust(bottom=0.2, top=0.85, left=0.1, right=0.95)

    os.makedirs("plots_output", exist_ok=True)
    save_path = f"plots_output/histogram_eta_{eta_str}_p_{p_target}.png"

    # Fixed size saving (removed bbox_inches='tight')
    plt.savefig(save_path, dpi=300, facecolor='white')
    print(f"\nSaved histogram to {save_path}")

    pass # plt.show disabled

# ==========================================================
# EXECUTION COMMAND
# ==========================================================
# plot_optimal_state_histogram(
    N_target=3,
    K_target=2,
    eta_vec_target=[0.6, 0.6],
    p_target=0.8,
    filename="optimized_states_common_eta.csv"
)

import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import pandas as pd
import numpy as np
import os
import re

def plot_final_research_map_refined(N_target, K_target, filename, gamma=2.5, save_path=None):
    if not os.path.exists(filename):
        print(f"Database file '{filename}' not found.")
        return

    # FORCE PYCHARM TO USE WHITE BACKGROUND
    plt.rcParams.update({
        'figure.facecolor': 'white',
        'axes.facecolor': 'white',
        'savefig.facecolor': 'white',
        'text.color': 'black',
        'axes.labelcolor': 'black',
        'xtick.color': 'black',
        'ytick.color': 'black',
        'axes.edgecolor': 'black'
    })

    df = pd.read_csv(filename)
    mask = (df['N'] == N_target) & (df['K'] == K_target)
    subset = df[mask].copy().dropna(subset=['coeffs'])
    if subset.empty: return

    def extract_first_float(s):
        clean_s = str(s).replace('np.float64(', '').replace(')', '')
        vals = re.findall(r"[-+]?(?:\d*\.\d+|\d+)(?:[eE][-+]?\d+)?", clean_s)
        return float(vals[0]) if vals else 0.0

    subset['eta_scalar'] = subset['eta_vec'].apply(extract_first_float)
    subset['p'] = pd.to_numeric(subset['p'], errors='coerce')

    def generate_combinations_local(K, N):
        from itertools import product
        return [comb for comb in product(range(N + 1), repeat=K) if sum(comb) <= N]

    combs = generate_combinations_local(K_target, N_target)
    ref_combs_w = [tuple([0]*K_target)]
    for i in range(K_target):
        mode_fock = [0]*K_target
        mode_fock[i] = N_target
        ref_combs_w.append(tuple(mode_fock))

    # --- HUMPHREYS ET AL. TRUE OPTIMAL STATE ---
    d = K_target
    alpha = 1.0 / np.sqrt(d + np.sqrt(d))
    beta = np.sqrt(1.0 - d * (alpha**2))

    ref_opt = np.zeros(len(combs))
    indices_opt = [i for i, c in enumerate(combs) if c in ref_combs_w]

    # The reference mode |0,0> gets the beta amplitude
    ref_opt[indices_opt[0]] = beta
    # The phase modes |N,0> and |0,N> get the alpha amplitude
    for idx in indices_opt[1:]:
        ref_opt[idx] = alpha

    coeff_matrix = np.array([np.array([float(x) for x in str(row['coeffs']).split(',')]) for _, row in subset.iterrows()])

    # Calculate fidelity against the TRUE optimal state
    fidelities_sq = (coeff_matrix @ ref_opt)**2

    fig, ax = plt.subplots(figsize=(12, 10))
    fig.patch.set_alpha(1.0)

    norm = mcolors.PowerNorm(gamma=gamma, vmin=0, vmax=1)

    sc = ax.scatter(subset['p'], subset['eta_scalar'], c=fidelities_sq, s=600,
                    cmap='RdYlBu_r', norm=norm, edgecolors='black', zorder=3)

    ax.set_ylabel(r"Transmission Parameter $\eta$", fontsize=15, labelpad=10)
    ax.set_xlabel(r"Mixture Parameter $p$", fontsize=15, labelpad=10)

    # SHIFTED Y-COORDINATES: 1.02 for Full Transmission, -0.02 for No Transmission
    ax.text(-0.03, 1.02, "Full Transmission", transform=ax.get_yaxis_transform(),
            ha='right', va='center', fontweight='bold', color='darkgreen', fontsize=12)
    ax.text(-0.03, -0.02, "No Transmission", transform=ax.get_yaxis_transform(),
            ha='right', va='center', fontweight='bold', color='darkred', fontsize=12)
    ax.text(0.0, -0.03, "Full Interaction", transform=ax.get_xaxis_transform(),
            ha='center', va='top', fontweight='bold', fontsize=12)
    ax.text(1.0, -0.03, "No Interaction", transform=ax.get_xaxis_transform(),
            ha='center', va='top', fontweight='bold', fontsize=12)

    # CLEAN RHS COLORBAR
    cbar = plt.colorbar(sc, pad=0.04)
    cbar.ax.set_title(r"$F = |\langle \psi_{opt} | \psi_{sim} \rangle|^2$", fontsize=16, pad=15)
    # Standard numerical ticks only, removing all the dense textual states
    cbar.set_ticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])

    # Removed the plt.title completely here!

    ax.set_xlim(-0.05, 1.05)
    ax.set_ylim(0.05, 1.05)

    plt.grid(True, linestyle=':', alpha=0.4)
    # Adjusted 'top' to 0.95 since the title is now gone
    plt.subplots_adjust(left=0.15, right=0.95, top=0.95, bottom=0.1)

    if save_path:
        plt.savefig(save_path, dpi=300, facecolor='white', transparent=False)
        print(f"  -> Saved heatmap to {save_path}")

    pass # plt.show disabled

# Execute Heatmap Plotting
output_dir = "plots_output"
os.makedirs(output_dir, exist_ok=True)

print("Plotting the full 2D Fidelity Heatmap...")
plot_final_research_map_refined(
    N_target=3, K_target=2, filename="optimized_states_common_eta.csv",
    save_path=f"{output_dir}/heatmap_phase_transition.png"
)

def optimize_initial_state_crb(N, K, eta, p):
    """
    Optimizes the initial probe state to MINIMIZE the Cramér-Rao Bound (CRB).
    CRB is computed as the Trace of the inverse of the QFIM.
    """
    D = N + 1
    combinations = generate_combinations(K, N)
    num_coeffs = len(combinations)

    # FIX: Ensure eta is a list/vector, as N_tot expects an eta for each mode
    eta_vec = eta if isinstance(eta, (list, tuple, np.ndarray)) else [eta] * K

    generators = []
    for i in range(K):
        ops = [qeye(D)] * K
        ops[i] = num(D)
        generators.append(tensor(ops))

    basis_kets = []
    for comb in combinations:
        state_modes = [basis(D, n) for n in comb]
        basis_kets.append(tensor(state_modes))

    def objective(coeffs):
        norm = np.linalg.norm(coeffs)
        if norm < 1e-10: return 1e6

        initial_ket = sum((coeffs[i] / norm) * basis_kets[i] for i in range(num_coeffs))

        # Pass eta_vec instead of eta
        rho_out = N_tot(initial_ket, eta_vec, p, N, K)
        vals, vecs = diagonalize(rho_out)
        qfim = calculate_QFIM(rho_out, vals, vecs, generators)

        # Add a tiny regularization to avoid singular matrix inversion errors
        reg_qfim = qfim + np.eye(K) * 1e-12
        inv_qfim = np.linalg.pinv(reg_qfim)

        # Goal: Minimize the trace of the inverse QFIM (the true multiparameter CRB)
        return np.trace(inv_qfim)

    start_coeffs = np.array([1.0 / np.sqrt(num_coeffs)] * num_coeffs)
    cons = ({'type': 'eq', 'fun': lambda x: np.sum(x**2) - 1})
    bnds = [(0, 1) for _ in range(num_coeffs)]

    res = minimize(objective, start_coeffs, method='SLSQP', bounds=bnds, constraints=cons)

    # We return the Trace(QFIM) of the optimal state to keep the database and plotting backwards-compatible
    opt_coeffs = res.x
    norm = np.linalg.norm(opt_coeffs)
    initial_ket = sum((opt_coeffs[i] / norm) * basis_kets[i] for i in range(num_coeffs))

    # Pass eta_vec instead of eta here as well
    rho_out = N_tot(initial_ket, eta_vec, p, N, K)
    vals, vecs = diagonalize(rho_out)
    opt_qfim = calculate_QFIM(rho_out, vals, vecs, generators)
    max_qfi = np.trace(opt_qfim)

    return res.x, combinations, max_qfi


import pandas as pd

def save_simulation_results(master_results, N, K, filename="metrology_database_crb.csv"):
    """
    Parses the master_results dictionary and saves it to a CSV file.
    Format of master_results: { eta_value: [ {'p': float, 'max_qfi': float, 'coeffs': list, 'combinations': list}, ... ] }
    """
    rows = []
    for eta, data_list in master_results.items():
        for data in data_list:
            # Convert coefficients list to a comma-separated string
            coeffs_str = ",".join(map(str, data['coeffs']))
            # Convert combinations list of tuples to a semicolon-separated string
            combs_str = ";".join([",".join(map(str, comb)) for comb in data['combinations']])

            rows.append({
                'N': N,
                'K': K,
                'eta': eta,
                'p': data['p'],
                'max_qfi': data['max_qfi'],
                'coeffs': coeffs_str,
                'combinations': combs_str
            })

    # Convert to DataFrame and save
    df = pd.DataFrame(rows)
    df.to_csv(filename, index=False)
    print(f"Results successfully saved to {filename}")

import pandas as pd
import numpy as np
import os
import csv

def run_and_save_simulation_crb(N, K, eta_values, p_values, filename="metrology_database_crb.csv"):
    # 1. Load existing progress to avoid re-running
    completed = set()
    if os.path.exists(filename):
        try:
            df_existing = pd.read_csv(filename)
            # Find all the existing runs for this N and K
            mask = (df_existing['N'] == N) & (df_existing['K'] == K)
            for _, row in df_existing[mask].iterrows():
                # Extract eta safely
                e_val = float(row['eta'] if 'eta' in df_existing.columns else row['eta_vec'])
                p_val = float(row['p'])
                # Store rounded values to avoid floating point mismatch
                completed.add((round(e_val, 3), round(p_val, 3)))
        except Exception as e:
            print(f"Warning: Could not read existing database. Error: {e}")

    # 2. Open file in APPEND mode ('a') so it never overwrites
    file_exists = os.path.isfile(filename)

    with open(filename, mode='a', newline='') as f:
        writer = csv.writer(f)

        # Write headers if it's a brand new file
        if not file_exists or os.path.getsize(filename) == 0:
            writer.writerow(['N', 'K', 'eta', 'p', 'max_qfi', 'coeffs', 'combinations'])

        # 3. Iterate and optimize
        for eta in eta_values:
            eta_display = round(eta, 2)
            print(f"\n>>> Processing eta = {eta_display} (CRB Optimization)...")

            for p in p_values:
                # CHECK IF ALREADY CALCULATED
                if (round(eta, 3), round(p, 3)) in completed:
                    print(f"    p={p:.2f}: Already exists in database. Skipping...")
                    continue

                print(f"    p={p:.2f}: Optimizing...")
                opt_coeffs, combinations, max_qfi = optimize_initial_state_crb(N, K, eta, p)

                # Format the lists into strings for the CSV
                coeffs_str = ",".join(map(str, opt_coeffs.tolist()))
                combs_str = ";".join([",".join(map(str, comb)) for comb in combinations])

                # Write the single row to the CSV instantly
                writer.writerow([N, K, eta, p, max_qfi, coeffs_str, combs_str])
                f.flush() # Force the operating system to save to disk IMMEDIATELY

                # Add to completed set so we don't accidentally run it twice in the same loop
                completed.add((round(eta, 3), round(p, 3)))

# Define parameters for the run
N_param, K_param = 3, 2
eta_list_crb = [0.8, 0.85, 0.9, 0.95, 1.0]

# Feed it the ENTIRE list. It will skip what you already have and calculate the missing gaps.
p_list_full = np.round(np.linspace(0.0, 1.0, 11), 1)

# Execute
run_and_save_simulation_crb(N_param, K_param, eta_list_crb, p_list_full)
print("\nDatabase update complete!")

import numpy as np
from qutip import *
from itertools import product

def compute_numerical_multipixel_limits(N, K):
    """
    Constructs the theoretical optimal noiseless state from Humphreys 2013
    and pushes it through the numerical QFIM functions to find the exact computed bounds.
    """
    D = N + 1

    # 1. Calculate theoretical amplitudes
    alpha = 1.0 / np.sqrt(K + np.sqrt(K))
    beta = 1.0 / np.sqrt(1 + np.sqrt(K))

    # 2. Build combinations and basis kets
    combinations = [c for c in product(range(N + 1), repeat=K) if sum(c) <= N]
    basis_kets = [tensor([basis(D, n) for n in comb]) for comb in combinations]

    # 3. Assign amplitudes to the correct basis states
    coeffs = np.zeros(len(combinations))

    # Reference mode has all N particles -> signal modes are (0,0)
    idx_ref = combinations.index(tuple([0] * K))
    coeffs[idx_ref] = beta

    # Signal modes have all N particles one by one
    for i in range(K):
        comb_i = tuple([N if j == i else 0 for j in range(K)])
        idx_sig = combinations.index(comb_i)
        coeffs[idx_sig] = alpha

    coeffs = coeffs / np.linalg.norm(coeffs) # Normalize
    initial_ket = sum(coeffs[i] * basis_kets[i] for i in range(len(combinations)))

    # 4. Set up generators and noiseless channel (eta = 1.0, p = 0.0)
    generators = []
    for i in range(K):
        ops = [qeye(D)] * K
        ops[i] = num(D)
        generators.append(tensor(ops))

    eta_vec = [1.0] * K
    rho_out = N_tot(initial_ket, eta_vec, 0.0, N, K)

    # 5. Compute QFIM using your pipeline
    vals, vecs = diagonalize(rho_out)
    qfim = calculate_QFIM(rho_out, vals, vecs, generators)

    multipixel_qfi = np.trace(qfim)
    multipixel_crb = (K**2) / multipixel_qfi

    return multipixel_qfi, multipixel_crb

# Compute and save as global variables for the plots
COMPUTED_MULTIPIXEL_QFI, COMPUTED_MULTIPIXEL_CRB = compute_numerical_multipixel_limits(3, 2)

# Calculate analytical for printing comparison
analytical_qfi = (4 * 3**2) / (2 * (1 + np.sqrt(2))**2)
analytical_crb = (2 * (1 + np.sqrt(2))**2) / (4 * 3**2)

print("--- Noiseless Limit Verification (eta=1.0) ---")
print(f"Computed Multipixel QFI   : {COMPUTED_MULTIPIXEL_QFI:.6f}")
print(f"Analytical Multipixel QFI : {analytical_qfi:.6f}")
print("----------------------------------------------")
print(f"Computed Multipixel CRB   : {COMPUTED_MULTIPIXEL_CRB:.6f}")
print(f"Analytical Multipixel CRB : {analytical_crb:.6f}")

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os
import re
from matplotlib.gridspec import GridSpec

def format_ket_latex(coeffs, combinations, threshold=1e-2):
    terms = []
    for c, comb in zip(coeffs, combinations):
        if np.abs(c) > threshold:
            basis_str = f"|{','.join(map(str, comb))}\\rangle"
            terms.append(f"{c:.2f}{basis_str}")

    if not terms: return "$0$"
    chunks = [terms[i:i + 5] for i in range(0, len(terms), 5)]
    lines = []
    for i, chunk in enumerate(chunks):
        if i == 0: lines.append(f"${' + '.join(chunk)}$")
        else: lines.append(f"$+ {' + '.join(chunk)}$")
    return "\n".join(lines)

def plot_from_database_crb(N_target, K_target, eta_vec_target, filename, save_path=None, forced_size=None, hide_legend=False, hide_benchmark_labels=False):
    if not os.path.exists(filename):
        return None

    plt.rcParams.update({'figure.facecolor': 'white', 'axes.facecolor': 'white', 'savefig.facecolor': 'white'})

    df = pd.read_csv(filename)
    target_scalar = float(eta_vec_target[0])

    def is_eta_match(s):
        clean_s = str(s).replace('np.float64(', '').replace(')', '')
        vals = re.findall(r"[-+]?(?:\d*\.\d+|\d+)(?:[eE][-+]?\d+)?", clean_s)
        if not vals: return False
        return np.isclose(float(vals[0]), target_scalar, atol=1e-3)

    p_targets = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]
    mask = (df['N'] == N_target) & (df['K'] == K_target) & \
           (df['eta'] if 'eta' in df.columns else df['eta_vec']).apply(is_eta_match) & \
           (df['p'].apply(lambda x: any(np.isclose(x, pt, atol=0.01) for pt in p_targets)))

    subset = df[mask].sort_values(by='p')
    if subset.empty: return None

    # --- IDENTIFY THE BEST OPTIMIZED STATE ---
    best_idx = subset['max_qfi'].idxmax()
    best_p = subset.loc[best_idx, 'p']
    best_qfi = subset.loc[best_idx, 'max_qfi']

    current_size = forced_size if forced_size else 8.0

    # WIDENED FIGURE: Added 3.5 to the width to act as a dedicated space for the RHS legend,
    # preventing the main plots from getting horizontally squished!
    fig = plt.figure(figsize=(current_size + 3.5, current_size))
    gs = GridSpec(2, 1, height_ratios=[4, 1], hspace=0.45)
    ax_qfi = fig.add_subplot(gs[0])
    ax_var = fig.add_subplot(gs[1])

    d, N = K_target, N_target
    sql_val = (d**2) / N
    noon_val = (d**3) / (N**2)

    # RESTORED: Exactly using your global notebook variables
    computed_multipixel_crb = COMPUTED_MULTIPIXEL_CRB
    computed_multipixel_qfi = COMPUTED_MULTIPIXEL_QFI

    computed_bound = sql_val
    for _, row in subset.iterrows():
        if np.isclose(row['p'], 0.0, atol=0.01) and row['max_qfi'] > 0:
            computed_bound = (d**2) / row['max_qfi']
            break

    best_crb = (d**2) / best_qfi

    # --- WIDEN MARGINS TO PUSH TEXT OUTWARDS ---
    lhs_bound = min(computed_multipixel_crb, computed_bound, best_crb)
    rhs_bound = max(sql_val, computed_bound, best_crb)

    margin = (rhs_bound - lhs_bound) * 0.35
    if margin == 0: margin = 0.1
    ax_var.set_xlim(lhs_bound - margin, rhs_bound + margin)
    ax_var.set_ylim(-1.0, 1.0)
    bar_height = 0.15
    ax_var.fill_between([lhs_bound, rhs_bound], -bar_height, bar_height, color='#E0E6ED', zorder=1)

    benchmarks = [computed_multipixel_crb, noon_val, sql_val]
    benchmark_labels = [
        rf"$\mathbf{{Multipixel: \approx {computed_multipixel_crb:.2f}}}$",
        rf"$\mathbf{{Heisenberg: d^3/N^2 \approx {noon_val:.2f}}}$",
        rf"$\mathbf{{SQL: d^2/N \approx {sql_val:.2f}}}$"
    ]

    colors = plt.cm.viridis(np.linspace(0, 1, len(subset)))
    qfis = subset['max_qfi'].values

    # --- DRAW THE BENCHMARKS ---
    for i, (val, label) in enumerate(zip(benchmarks, benchmark_labels)):
        ax_var.plot([val, val], [-bar_height, bar_height], color='black', linewidth=2.5, zorder=2)
        if not hide_benchmark_labels:
            if i == 0:
                ax_var.text(val - margin*0.02, bar_height + 0.1, label, color='black', fontsize=11, ha='right', va='bottom')
            elif i == 2:
                ax_var.text(val + margin*0.02, bar_height + 0.1, label, color='black', fontsize=11, ha='left', va='bottom')
            else:
                ax_var.text(val, bar_height + 0.1, label, color='black', fontsize=11, ha='center', va='bottom')

    # --- PLOT SCATTER ---
    for i, (idx, row) in enumerate(subset.iterrows()):
        p_val = row['p']
        qfi_val = row['max_qfi']
        var_val = (d**2) / qfi_val
        coeffs = [float(x) for x in str(row['coeffs']).split(',')]
        combs = [tuple(map(int, c.split(','))) for c in str(row['combinations']).split(';')]

        latex_ket = format_ket_latex(coeffs, combs)
        legend_label = f"p={p_val:.2f}: {latex_ket}"

        ax_qfi.scatter(p_val, qfi_val, label=legend_label, s=120, color=colors[i], edgecolors='black', zorder=5)

        if np.isclose(p_val, best_p, atol=0.01):
            ax_var.plot([var_val, var_val], [-bar_height, bar_height], color=colors[i], linewidth=8, solid_capstyle='butt', zorder=5)
            ax_var.text(var_val, -bar_height - 0.15, rf"Best Opt: {var_val:.2f}", color=colors[i], fontsize=11, ha='center', va='top', weight='bold')

    hl_qfi_val = (d**2) / noon_val

    ax_qfi.axhline(y=hl_qfi_val, color='red', linestyle='--', linewidth=2, zorder=2, label=rf"Heisenberg Limit QFI ($\approx {hl_qfi_val:.2f}$)")
    ax_qfi.axhline(y=computed_multipixel_qfi, color='dodgerblue', linestyle='-.', linewidth=2, zorder=2, label=rf"Multipixel QFI ($\approx {computed_multipixel_qfi:.2f}$)")

    # --- EXTEND Y-AXIS DRAMATICALLY FOR LEGEND CLEARANCE ---
    all_qfi_values = list(qfis) + [hl_qfi_val, computed_multipixel_qfi]
    if len(all_qfi_values) > 0:
        q_min, q_max = np.min(all_qfi_values), np.max(all_qfi_values)
        if q_max == q_min: q_min -= 1; q_max += 1
        ax_qfi.set_ylim(q_min - 0.1*(q_max-q_min), q_max + 0.5*(q_max-q_min))

    ax_qfi.set_xlabel(r"Mixture Parameter $p$ (0 = Full Interaction, 1 = No Interaction)", fontsize=11)
    ax_qfi.set_ylabel(r"Maximized $\mathrm{Tr}(\mathcal{I}_{\boldsymbol{\theta}})$", fontsize=11)
    ax_qfi.grid(True, linestyle=':', alpha=0.4)

    ax_var.set_yticks([]); ax_var.spines['left'].set_visible(False)
    ax_var.spines['right'].set_visible(False); ax_var.spines['top'].set_visible(False)
    ax_var.spines['bottom'].set_visible(False); ax_var.set_xticks([])
    ax_var.set_title(rf"Cramer-Rao Bound of State at $p={best_p:.2f}$", fontsize=13, pad=10)

    # Limit legends (Top Left inside the plot)
    handles_qfi, labels_qfi = ax_qfi.get_legend_handles_labels()
    line_handles = [h for h, l in zip(handles_qfi, labels_qfi) if "QFI" in l]
    line_labels = [l for l in labels_qfi if "QFI" in l]

    limit_legend = ax_qfi.legend(line_handles, line_labels, loc='upper left', fontsize=9, frameon=True, facecolor='white', framealpha=0.9)
    ax_qfi.add_artist(limit_legend)

    # --- FIXED RIGHT-HAND SIDE LEGEND ---
    if not hide_legend:
        handles, labels = ax_qfi.get_legend_handles_labels()
        state_handles_labels = [(h, l) for h, l in zip(handles, labels) if "QFI" not in l]
        if state_handles_labels:
            sh, sl = zip(*state_handles_labels)

            # Shrunk the right bound to exactly 0.60 so the plot leaves space for the legend.
            # Kept top=0.95, which is exactly where we will pin the legend.
            plt.subplots_adjust(left=0.10, bottom=0.15, right=0.60, top=0.95)

            # Switched to fig.legend so it never gets buried or overridden by the axis.
            # bbox_to_anchor=(0.62, 0.95) places it right next to the plot border (0.60)
            # and perfectly aligns the top edges (0.95).
            states_legend = fig.legend(sh, sl, title=rf"Optimized States $|\psi\rangle$",
                                       loc='upper left', bbox_to_anchor=(0.62, 0.95),
                                       fontsize=10, frameon=True, facecolor='white',
                                       ncol=1, title_fontsize='12', labelspacing=2.0)
    else:
        plt.subplots_adjust(left=0.15, bottom=0.15, right=0.95, top=0.95)

    if forced_size is None:
        fig.canvas.draw()
        bbox = fig.get_tightbbox(fig.canvas.get_renderer())
        plt.close(fig)
        return bbox.width, bbox.height
    else:
        if save_path:
            plt.savefig(save_path, dpi=300, facecolor='white', transparent=False, bbox_inches='tight')
            print(f"  -> Saved figure to {save_path}")
        pass # plt.show disabled
        return None

# Execution Block for 1D Plots
