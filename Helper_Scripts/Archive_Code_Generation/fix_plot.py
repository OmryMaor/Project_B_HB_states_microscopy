import json

notebook_path = "c:/Users/omrym/Technion/Semester9/Project B/HB_states_ProjectB/HB_QFI_loss.ipynb"

with open(notebook_path, 'r', encoding='utf-8') as f:
    nb = json.load(f)

for cell in nb.get('cells', []):
    if cell['cell_type'] == 'code':
        # Modify the all etas cell properly
        new_source = []
        for line in cell['source']:
            if "linewidth=2" in line and "linewidth=2, color=colors" in line:
                line = line.replace("linewidth=2", "linestyle='None'")
            new_source.append(line)
        cell['source'] = new_source

with open(notebook_path, 'w', encoding='utf-8') as f:
    json.dump(nb, f, indent=1)

print("Notebook plot fixed.")
