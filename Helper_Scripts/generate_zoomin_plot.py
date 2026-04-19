import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import os

# ---- Output folder ----
out_folder = r"c:\Users\omrym\Technion\Semester9\Project B\HB_states_ProjectB\plots_output\Humphreys_Replications"
os.makedirs(out_folder, exist_ok=True)

# Data collected from previous run (d=2)
N_values_full = [3, 6, 9, 12, 15, 18]
N_values_zoom = [9, 12, 15, 18]

# Pure CRB (eta=1.0)
crbs_pure = [0.5000, 0.1667, 0.0833, 0.0500, 0.0333, 0.0238]
crbs_pure_zoom = crbs_pure[2:]

# Mixed CRB
mixed_data = {
    0.95: [0.5270, 0.1803, 0.0926, 0.0570, 0.0390, 0.0286],
    0.90: [0.5585, 0.1966, 0.1038, 0.0658, 0.0463, 0.0348],
    0.85: [0.5956, 0.2161, 0.1176, 0.0766, 0.0553, 0.0426]
}

etas = [0.95, 0.90, 0.85]
colors = ['purple', 'orange', 'brown']

# Benchmarks for d=2
d_fixed = 2
multipixel_zoom = [(d_fixed * (1 + np.sqrt(d_fixed))**2) / (4 * N**2) for N in N_values_zoom]
noon_zoom       = [(d_fixed**3) / N**2 for N in N_values_zoom]

# ---- Zoom-in Plot ----
fig, ax = plt.subplots(figsize=(8, 5))

# Pure state baseline
ax.plot(N_values_zoom, crbs_pure_zoom, 's-', color='green', markersize=8,
        label=r'HB(n,%d) $\eta=1.0$' % (d_fixed + 1))

# Mixed state curves
for i, eta in enumerate(etas):
    ax.plot(N_values_zoom, mixed_data[eta][2:], 'o-', color=colors[i], markersize=8,
            label=r'HB(n,%d) $\eta=%s$' % (d_fixed + 1, eta))

# Ideal limits
ax.plot(N_values_zoom, multipixel_zoom, '--', color='red', lw=2,
        label=r'Optimal $|\psi_s\rangle$ (Ideal)')
ax.plot(N_values_zoom, noon_zoom, '--', color='blue', lw=2,
        label='N00N States (Ideal)')

ax.set_xlabel(f"Total photon number N = {d_fixed+1}n", fontsize=12)
ax.set_ylabel(r"Total variance $|\Delta\theta|^2$", fontsize=12)
ax.set_title(f"CRB vs N for d={d_fixed} (Zoom-in)", fontsize=14)
ax.set_xticks(N_values_zoom)
ax.set_ylim(0, 0.2)
ax.grid(True, ls=':', alpha=0.6)
ax.legend(fontsize=10)

path = os.path.join(out_folder, "QCRB_vs_N_d2_zoomin.png")
fig.savefig(path, dpi=300, bbox_inches='tight')
plt.close(fig)
print(f"Saved -> {path}")
