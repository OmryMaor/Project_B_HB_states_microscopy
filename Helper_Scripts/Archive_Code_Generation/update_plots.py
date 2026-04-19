import json

notebook_path = "c:/Users/omrym/Technion/Semester9/Project B/HB_states_ProjectB/HB_QFI_loss.ipynb"

with open(notebook_path, 'r', encoding='utf-8') as f:
    nb = json.load(f)

for cell in nb.get('cells', []):
    if cell['cell_type'] == 'code':
        # Modify the all etas cell (cell -3)
        if "def plot_hb_all_etas" in "".join(cell['source']):
            new_source = []
            for line in cell['source']:
                if "ax_qfi.plot(subset_hb['p']" in line and "linewidth=2" in line:
                    line = line.replace("linewidth=2", "linestyle='None'")
                new_source.append(line)
            cell['source'] = new_source
            
        # Modify the N dependence cell (cell -2)
        if "CRB vs N" in "".join(cell['source']):
            new_source = []
            for line in cell['source']:
                if "plt.plot(N_values, crb_N_results" in line and "linewidth=2" in line:
                    line = line.replace("linewidth=2", "linestyle='None'")
                if "CRB_vs_N_d{d_fixed}_eta{eta_target}_p{p_target}.png" in line:
                    line = line.replace(".png", "_scatter.png")
                new_source.append(line)
            cell['source'] = new_source
            
        # Modify the d dependence cell (cell -1)
        if "CRB vs d" in "".join(cell['source']):
            new_source = []
            for line in cell['source']:
                if "d_values = [2, 3]" in line:
                    line = "d_values = [1, 2, 3, 4]\n"
                if "plt.plot(d_values, crb_d_results" in line and "linewidth=2" in line:
                    line = line.replace("linewidth=2", "linestyle='None'")
                if "CRB_vs_d_N{N_fixed}_eta{eta_target}_p{p_target}.png" in line:
                    line = line.replace(".png", "_scatter.png")
                new_source.append(line)
            cell['source'] = new_source

with open(notebook_path, 'w', encoding='utf-8') as f:
    json.dump(nb, f, indent=1)

print("Notebook updated successfully.")
