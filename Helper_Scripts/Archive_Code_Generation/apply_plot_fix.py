import re

with open('create_notebook.py', 'r', encoding='utf-8') as f:
    text = f.read()

prefix = "evals = r'''"
suffix = "'''\nnb.cells.append(new_code_cell(evals))"

start_idx = text.find(prefix)
end_idx = text.find(suffix, start_idx) + len(suffix)

new_evals = r"""
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
            terms.append(f"{c:.2f}|{','.join(map(str, comb))}\\rangle")
    return " + ".join(terms)

def plot_hb_vs_optimized(df_hb, df_opt, N_target, K_target, eta_target, output_dir):
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
                latex_ket = f"$p={row['p']:.1f}: $ {format_ket_latex(c_arr, cb_arr)}"
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
    
    if ket_labels_handles:
        l_hand = [h[0] for h in ket_labels_handles]
        l_lab = [h[1] for h in ket_labels_handles]
        
        hb_ket_str = r"$|HB\rangle: $ Fixed input superposition" 
        l_hand.append(hb_h)
        l_lab.append(hb_ket_str)
        
        leg2 = plt.legend(l_hand, l_lab, loc='center left', bbox_to_anchor=(1.02, 0.5), title=r"$\mathbf{Optimized\ States\ (\eta=" + str(eta_target) + ")}$",
                 fontsize=9, frameon=True)
    
    plt.subplots_adjust(left=0.08, bottom=0.12, right=0.6, top=0.92)
    
    import os
    os.makedirs(output_dir, exist_ok=True)
    filename = f"{output_dir}/HB_vs_Optimized_eta_{eta_target}.png"
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    print(f"Saved aesthetically matched scatter plot to {filename}")
    plt.close()
"""

replacement = prefix + "\n" + new_evals + "\n" + suffix
new_text = text[:start_idx] + replacement + text[end_idx:]

with open('create_notebook.py', 'w', encoding='utf-8') as f:
    f.write(new_text)

print('Success')
