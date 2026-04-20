import numpy as np
import matplotlib.pyplot as plt
import os

# --- Setup ---
out_folder = "plots_output/Humphreys_Replications"
os.makedirs(out_folder, exist_ok=True)

# --- Data from the 30-minute run ---
N_vals_a = [3, 6, 9, 12, 15, 18]
hb_a = {
    1.0: [0.5000, 0.1667, 0.0833, 0.0500, 0.0333, 0.0238],
    0.95: [0.5270, 0.1803, 0.0926, 0.0570, 0.0390, 0.0286],
    0.90: [0.5585, 0.1966, 0.1038, 0.0658, 0.0463, 0.0348],
    0.85: [0.5956, 0.2161, 0.1176, 0.0766, 0.0553, 0.0426]
}

d_vals_b = [1, 2, 3, 4]
hb_b = {
    1.0: [0.2500, 0.5000, 0.7500, 1.0000],
    0.95: [0.2635, 0.5270, 0.7907, 1.0543],
    0.90: [0.2793, 0.5585, 0.8385, 1.1184],
    0.85: [0.2980, 0.5956, 0.8951, 1.1946]
}

etas = [0.95, 0.9, 0.85]
colors = ['purple', 'orange', 'brown']

def noisy_noon(N, d, eta):
    n_mode = N / d
    return (d**3 * (1 + eta**n_mode)) / (2 * (N**2) * (eta**n_mode))

# --- PLOT A: QCRB vs N ---
fig, ax = plt.subplots(figsize=(8, 5))
ax.plot(N_vals_a, hb_a[1.0], 's-', color='green', markersize=8, label=r'HB $\eta=1.0$')
for i, eta in enumerate(etas):
    ax.plot(N_vals_a, hb_a[eta], 'o-', color=colors[i], markersize=8, label=r'HB $\eta=%s$' % eta)
    noisy_n = [noisy_noon(N, 2, eta) for N in N_vals_a]
    ax.plot(N_vals_a, noisy_n, '--', color=colors[i], alpha=0.7, label=r'N00N $\eta=%s$' % eta)

# Ideal benchmarks
multipixel_a = [(2 * (1 + np.sqrt(2))**2) / (4 * N**2) for N in N_vals_a]
noon_ideal_a = [(2**3) / N**2 for N in N_vals_a]
ax.plot(N_vals_a, multipixel_a, ':', color='red', lw=1.5, label='Optimal (Ideal)')
ax.plot(N_vals_a, noon_ideal_a, '--', color='green', lw=1.5, label='N00N (Ideal)')

ax.set_xlabel('Total photon number N = 3n')
ax.set_ylabel(r'Total variance $|\Delta\theta|^2$')
ax.set_title('QCRB vs N (d=2)')
ax.set_xticks(N_vals_a)
ax.grid(True, ls=':', alpha=0.6)
ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()
fig.savefig(f"{out_folder}/QCRB_vs_N_d2_noise.png", dpi=300, bbox_inches='tight')

# --- PLOT A ZOOM-IN ---
ax.set_ylim(0, 0.2)
ax.set_xlim(8, 19)
ax.set_title('QCRB vs N (d=2) - Zoom-in')
fig.savefig(f"{out_folder}/QCRB_vs_N_d2_zoomin.png", dpi=300, bbox_inches='tight')
plt.close(fig)

# --- PLOT B: QCRB vs d ---
fig, ax = plt.subplots(figsize=(8, 5))
ax.plot(d_vals_b, hb_b[1.0], 's-', color='green', markersize=8, label=r'HB $\eta=1.0$')
for i, eta in enumerate(etas):
    ax.plot(d_vals_b, hb_b[eta], 'o-', color=colors[i], markersize=8, label=r'HB $\eta=%s$' % eta)
    # N = d+1
    noisy_n = [noisy_noon(d+1, d, eta) for d in d_vals_b]
    ax.plot(d_vals_b, noisy_n, '--', color=colors[i], alpha=0.7, label=r'N00N $\eta=%s$' % eta)

multipixel_b = [(d * (1 + np.sqrt(d))**2) / (4 * (d+1)**2) for d in d_vals_b]
noon_ideal_b = [(d**3) / (d+1)**2 for d in d_vals_b]
ax.plot(d_vals_b, multipixel_b, ':', color='red', lw=1.5, label='Optimal (Ideal)')
ax.plot(d_vals_b, noon_ideal_b, '--', color='green', lw=1.5, label='N00N (Ideal)')

ax.set_xlabel('Number of phases d (N = d+1)')
ax.set_ylabel(r'Total variance $|\Delta\theta|^2$')
ax.set_title('QCRB vs d (n=1)')
ax.set_xticks(d_vals_b)
ax.grid(True, ls=':', alpha=0.6)
ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()
fig.savefig(f"{out_folder}/QCRB_vs_d_noise.png", dpi=300, bbox_inches='tight')
plt.close(fig)

print("Regenerated all plots with corrected colors and zoom-in benchmarks.")
