import json
import os

notebook_path = "c:/Users/omrym/Technion/Semester9/Project B/HB_states_ProjectB/HB_QFI_loss.ipynb"

with open(notebook_path, 'r', encoding='utf-8') as f:
    nb = json.load(f)

new_cells = []

# Cell 1: Function to plot all etas together
code1 = """
# == NEW PLOT: All etas in one plot ==

def plot_hb_all_etas(df_hb, df_opt, hb_coeffs, combs, N_target, K_target, eta_values, output_dir, file_prefix=\"HB_vs_Optimized_All_Etas\"):
    plt.rcParams.update({'figure.facecolor': 'white', 'axes.facecolor': 'white', 'savefig.facecolor': 'white'})
    
    fig, ax_qfi = plt.subplots(figsize=(10.5, 6.5))
    
    colors = plt.cm.plasma(np.linspace(0, 0.9, len(eta_values)))
    
    for idx, eta in enumerate(eta_values):
        subset_hb = df_hb[np.isclose(df_hb['eta'], eta, atol=1e-3)].sort_values(by='p')
        if not subset_hb.empty:
            ax_qfi.plot(subset_hb['p'], subset_hb['total_var'], marker='o', markersize=8, 
                        linewidth=2, color=colors[idx], label=f"HB $\\\\eta={eta}$")
            
        if df_opt is not None:
            subset_opt = df_opt[np.isclose(df_opt['eta'], eta, atol=1e-3)].sort_values(by='p')
            valid_ps = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]
            subset_opt = subset_opt[pd.to_numeric(subset_opt['p']).apply(lambda x: any(np.isclose(x, v, atol=1e-3) for v in valid_ps))]
            if not subset_opt.empty:
                ax_qfi.scatter(subset_opt['p'], subset_opt['total_var'], marker='*', s=200, 
                               edgecolors='black', color=colors[idx], alpha=0.9, label=f"Opt $\\\\eta={eta}$")

    d, N = K_target, N_target
    noon_val = (d**3) / (N**2)
    computed_multipixel_crb = (d * (1 + np.sqrt(d))**2) / (4 * (N**2))
    
    ax_qfi.axhline(y=noon_val, color='red', linestyle='--', linewidth=2, zorder=1, label=f"Heisenberg Limit ($\\\\approx {noon_val:.2f}$)")
    ax_qfi.axhline(y=computed_multipixel_crb, color='dodgerblue', linestyle='-.', linewidth=2, zorder=1, label=f"Multipixel CRB ($\\\\approx {computed_multipixel_crb:.2f}$)")

    ax_qfi.set_xlabel(\"Mixture Parameter $p$ (0 = Full Interaction, 1 = No Interaction)\", fontsize=11)
    ax_qfi.set_ylabel(r\"Minimized $\\mathrm{Tr}(\\mathcal{I}_{\\boldsymbol{\\theta}}^{-1})$\", fontsize=11)
    ax_qfi.grid(True, linestyle=':', alpha=0.4)
    ax_qfi.set_title(f\"State Comparison for N={N}, d={d} (All $\\eta$)\", fontsize=13, pad=10)
    
    handles, labels = ax_qfi.get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    ax_qfi.legend(by_label.values(), by_label.keys(), loc='best', fontsize=10, frameon=True, facecolor='white', bbox_to_anchor=(1.05, 1))
    
    plt.subplots_adjust(left=0.08, bottom=0.12, right=0.75, top=0.92)
    
    os.makedirs(output_dir, exist_ok=True)
    filename = f\"{output_dir}/{file_prefix}_N{N_target}_d{K_target}.png\"
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    print(f\"Saved all-etas plot to {filename}\")
    pass # plt.show disabled

# Generate the plot for the initial case N=3, d=2
plot_hb_all_etas(df_hb, df_opt, hb_coeffs, combs, N_photons, d_phases, eta_grid, out_folder)

# And for higher dimensions
for (N_target, d_target) in configs:
    # We re-calculate since df_hb_h was overwritten in the loop. 
    # Or actually wait, we can just run it
    print(f\"Recalculating grid for N={N_target}, d={d_target} to plot all etas...\")
    hb_c, b_kets, cmbs = initialize_HB_state(N_target, d_target)
    df_hb_h = evaluate_grid(N_target, d_target, hb_c, b_kets, eta_grid, p_grid)
    plot_hb_all_etas(df_hb_h, None, hb_c, cmbs, N_target, d_target, eta_grid, out_folder_high, file_prefix=\"HB_alone_All_Etas\")
"""

new_cells.append({
    "cell_type": "code",
    "execution_count": None,
    "metadata": {},
    "outputs": [],
    "source": [line + "\n" for line in code1.split('\n')[:-1]] + [code1.split('\n')[-1]]
})


# Cell 2: N dependence for eta=0.95, p=1.0
code2 = """
# == NEW PLOT: CRB vs N (number of photons) for fixed d, eta=0.95, p=1.0 ==

d_fixed = 3
eta_target = 0.95
p_target = 1.0

N_values = [2, 3, 4, 5, 6, 7]
crb_N_results = []
qfi_N_results = []

print(f\"Evaluating N-dependence for d={d_fixed}, eta={eta_target}, p={p_target}...\")
for N_test in N_values:
    # Initialize state
    hb_c, b_kets, cmbs = initialize_HB_state(N_test, d_fixed)
    initial_ket = construct_ket(hb_c, b_kets)
    
    # Evaluate metric
    qfi, crb = evaluate_state_metrics(initial_ket, N_test, d_fixed, [eta_target]*d_fixed, p_target)
    crb_N_results.append(crb)
    qfi_N_results.append(qfi)
    print(f\"N={N_test} | QFI={qfi:.3f}, CRB={crb:.3f}\")

plt.figure(figsize=(8, 5))
plt.plot(N_values, crb_N_results, marker='o', linewidth=2, color='darkorange', label=f'HB CRB ($\\eta={eta_target}$, p={p_target})')

# Also plot the multipixel CRB and HL for reference
multipixel_crbs = [(d_fixed * (1 + np.sqrt(d_fixed))**2) / (4 * (n**2)) for n in N_values]
heisenberg_limits = [(d_fixed**3) / (n**2) for n in N_values]

plt.plot(N_values, multipixel_crbs, linestyle='-.', color='dodgerblue', label='Multipixel CRB (Ideal)')
plt.plot(N_values, heisenberg_limits, linestyle='--', color='red', label='Heisenberg Limit')

plt.xlabel(\"Number of Photons (N)\", fontsize=12)
plt.ylabel(r\"Minimum Variance (CRB)\", fontsize=12)
plt.title(f\"CRB vs N for d={d_fixed} ($\\\\eta={eta_target}, p={p_target}$)\", fontsize=14)
plt.grid(True, linestyle=':', alpha=0.6)
plt.legend()

filename_N = f\"{out_folder_high}/CRB_vs_N_d{d_fixed}_eta{eta_target}_p{p_target}.png\"
plt.savefig(filename_N, dpi=300, bbox_inches='tight')
print(f\"Saved N-dependence plot to {filename_N}\")
pass # plt.show disabled
"""
new_cells.append({
    "cell_type": "code",
    "execution_count": None,
    "metadata": {},
    "outputs": [],
    "source": [line + "\n" for line in code2.split('\n')[:-1]] + [code2.split('\n')[-1]]
})

# Cell 3: d dependence for eta=0.95, p=1.0
code3 = """
# == NEW PLOT: CRB vs d (number of pixels/phases) for fixed N, eta=0.95, p=1.0 ==

N_fixed = 6
eta_target = 0.95
p_target = 1.0

d_values = [2, 3, 4]  # d=1 isn't multi-phase; d=5 might take too long, so let's stick to 2,3,4. 
crb_d_results = []

print(f\"Evaluating d-dependence for N={N_fixed}, eta={eta_target}, p={p_target}...\")
for d_test in d_values:
    # Initialize state
    hb_c, b_kets, cmbs = initialize_HB_state(N_fixed, d_test)
    initial_ket = construct_ket(hb_c, b_kets)
    
    # Evaluate metric
    qfi, crb = evaluate_state_metrics(initial_ket, N_fixed, d_test, [eta_target]*d_test, p_target)
    crb_d_results.append(crb)
    print(f\"d={d_test} | QFI={qfi:.3f}, CRB={crb:.3f}\")

plt.figure(figsize=(8, 5))
plt.plot(d_values, crb_d_results, marker='s', linewidth=2, color='purple', label=f'HB CRB ($\\eta={eta_target}$, p={p_target})')

multipixel_crbs_d = [(d * (1 + np.sqrt(d))**2) / (4 * (N_fixed**2)) for d in d_values]
heisenberg_limits_d = [(d**3) / (N_fixed**2) for d in d_values]

plt.plot(d_values, multipixel_crbs_d, linestyle='-.', color='dodgerblue', label='Multipixel CRB (Ideal)')
plt.plot(d_values, heisenberg_limits_d, linestyle='--', color='red', label='Heisenberg Limit')

plt.xlabel(\"Number of Pixels/Phases (d)\", fontsize=12)
plt.ylabel(r\"Minimum Variance (CRB)\", fontsize=12)
plt.title(f\"CRB vs d for N={N_fixed} ($\\\\eta={eta_target}, p={p_target}$)\", fontsize=14)
plt.xticks(d_values)
plt.grid(True, linestyle=':', alpha=0.6)
plt.legend()

filename_d = f\"{out_folder_high}/CRB_vs_d_N{N_fixed}_eta{eta_target}_p{p_target}.png\"
plt.savefig(filename_d, dpi=300, bbox_inches='tight')
print(f\"Saved d-dependence plot to {filename_d}\")
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

print("Cells successfully appended to notebook!")
