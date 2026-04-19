import json

notebook_path = "c:/Users/omrym/Technion/Semester9/Project B/HB_states_ProjectB/HB_QFI_loss.ipynb"

with open(notebook_path, 'r', encoding='utf-8') as f:
    nb = json.load(f)

for cell in nb.get('cells', []):
    if cell['cell_type'] == 'code':
        # Replace N_tot
        for i, line in enumerate(cell['source']):
            if "def N_tot(" in line:
                # The function is the last one in this cell
                pass
            if "return p * N_free(" in line and "(1.0 - p) * N_int(initial_ket" in line:
                cell['source'][i] = "    if p == 1.0:\n        return N_free(initial_ket, eta_vec, N, K)\n    elif p == 0.0:\n        return N_int(initial_ket, eta_vec, N, K)\n    return p * N_free(initial_ket, eta_vec, N, K) + (1.0 - p) * N_int(initial_ket, eta_vec, N, K)\n"
        
        # Replace d_values = [2, 3, 4] with d_values = [2, 3] to be safe
        for i, line in enumerate(cell['source']):
            if "d_values = [2, 3, 4]" in line:
                cell['source'][i] = line.replace("d_values = [2, 3, 4]", "d_values = [2, 3]")

with open(notebook_path, 'w', encoding='utf-8') as f:
    json.dump(nb, f, indent=1)

print("Notebook fixes applied.")
