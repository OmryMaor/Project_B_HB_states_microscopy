import json
notebook_path = 'c:/Users/omrym/Technion/Semester9/Project B/HB_states_ProjectB/HB_QFI_loss.ipynb'
with open(notebook_path, 'r', encoding='utf-8') as f: nb = json.load(f)
for cell in nb['cells']:
    if cell['cell_type'] == 'code':
        new_source = []
        for line in cell['source']:
            line = line.replace('f"{out_folder_high}/Replication_Fig', 'f"plots_output/Humphreys_Replications/Replication_Fig')
            new_source.append(line)
        cell['source'] = new_source
with open(notebook_path, 'w', encoding='utf-8') as f: json.dump(nb, f, indent=1)
