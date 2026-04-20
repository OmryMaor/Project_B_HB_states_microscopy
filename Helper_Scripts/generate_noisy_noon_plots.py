import numpy as np
import matplotlib.pyplot as plt
import os

# --- Setup Paths ---
out_folder = "plots_output/Humphreys_Replications"
os.makedirs(out_folder, exist_ok=True)

# --- Data from compute_N12.py (Latest results) ---
etas = [0.95, 0.90, 0.85]
colors_eta = ['purple', 'orange', 'brown']

# Fig 3A: d=3 fixed, n=1,2,3
d_a = 3
N_val_a = np.array([4, 8, 12])
# HB QCRB values (including the N=12 computed yesterday)
hb_a = {
    0.95: [0.7907, 0.2706, 0.1384],
    0.90: [0.8385, 0.2953, 0.1557],
    0.85: [0.8951, 0.3252, 0.1782],
}
hb_pure_a = [0.7500, 0.2500, 0.1250]

# Fig 3B: n=1 fixed, d=1,2,3,4
d_val_b = np.array([1, 2, 3, 4])
N_val_b = d_val_b + 1
hb_b = {
    0.95: [0.2635, 0.5270, 0.7907, 1.0543],
    0.90: [0.2793, 0.5585, 0.8385, 1.1184],
    0.85: [0.2980, 0.5956, 0.8951, 1.1946],
}
hb_pure_b = [0.2500, 0.5000, 0.7500, 1.0000]

# --- Benchmarking Functions ---
def noisy_noon(N, d, eta):
    # Formula: d^3 * (1 + eta^(N/d)) / (2 * N^2 * eta^(N/d))
    # N/d is photons per mode
    n_mode = N / d
    return (d**3 * (1 + eta**n_mode)) / (2 * (N**2) * (eta**n_mode))

def ideal_multipixel(N, d):
    return (d * (1 + np.sqrt(d))**2) / (4 * (N**2))

def ideal_noon(N, d):
    return (d**3) / (N**2)

# ===================================================================
# PLOT A: QCRB vs N (d=3)
# ===================================================================
plt.figure(figsize=(9, 6))

# 1. Pure HB
plt.plot(N_val_a, hb_pure_a, 's-', color='green', markersize=8, label=r'HB(n, 4) $\eta=1.0$')

# 2. Noisy curves
for i, eta in enumerate(etas):
    # HB State (Solid)
    plt.plot(N_val_a, hb_a[eta], 'o-', color=colors_eta[i], markersize=8, 
             label=r'HB(n, 4) $\eta=%s$' % eta)
    
    # Noisy N00N (Dashed, same color)
    noisy_n = [noisy_noon(N, d_a, eta) for N in N_val_a]
    plt.plot(N_val_a, noisy_n, '--', color=colors_eta[i], alpha=0.7, 
             label=r'N00N $\eta=%s$' % eta)

# 3. Ideal benchmarks (Neutral)
plt.plot(N_val_a, [ideal_multipixel(N, d_a) for N in N_val_a], ':', color='red', lw=1.5, label='Optimal (Ideal)')
plt.plot(N_val_a, [ideal_noon(N, d_a) for N in N_val_a], ':', color='blue', lw=1.5, label='N00N (Ideal)')

plt.xlabel('Total Photon Number N (4n)', fontsize=12)
plt.ylabel(r'Total Variance $|\Delta\theta|^2$', fontsize=12)
plt.title('HB vs Noisy N00N: QCRB vs N (d=3)', fontsize=14)
plt.xticks(N_val_a)
plt.grid(True, ls=':', alpha=0.6)
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=9)
plt.tight_layout()
plt.savefig(f"{out_folder}/QCRB_vs_N_noisy_noon_comparison.png", dpi=300)
plt.close()

# ===================================================================
# PLOT B: QCRB vs d (n=1)
# ===================================================================
plt.figure(figsize=(9, 6))

# 1. Pure HB
plt.plot(d_val_b, hb_pure_b, 's-', color='green', markersize=8, label=r'HB(1, d) $\eta=1.0$')

# 2. Noisy curves
for i, eta in enumerate(etas):
    # HB State (Solid)
    plt.plot(d_val_b, hb_b[eta], 'o-', color=colors_eta[i], markersize=8, 
             label=r'HB(1, d) $\eta=%s$' % eta)
    
    # Noisy N00N (Dashed, same color)
    noisy_n = [noisy_noon(N, d, eta) for N, d in zip(N_val_b, d_val_b)]
    plt.plot(d_val_b, noisy_n, '--', color=colors_eta[i], alpha=0.7, 
             label=r'N00N $\eta=%s$' % eta)

# 3. Ideal benchmarks (Neutral)
plt.plot(d_val_b, [ideal_multipixel(N, d) for N, d in zip(N_val_b, d_val_b)], ':', color='red', lw=1.5, label='Optimal (Ideal)')
plt.plot(d_val_b, [ideal_noon(N, d) for N, d in zip(N_val_b, d_val_b)], ':', color='blue', lw=1.5, label='N00N (Ideal)')

plt.xlabel('Number of Phases d (N = d+1)', fontsize=12)
plt.ylabel(r'Total Variance $|\Delta\theta|^2$', fontsize=12)
plt.title('HB vs Noisy N00N: QCRB vs d (n=1)', fontsize=14)
plt.xticks(d_val_b)
plt.grid(True, ls=':', alpha=0.6)
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=9)
plt.tight_layout()
plt.savefig(f"{out_folder}/QCRB_vs_d_noisy_noon_comparison.png", dpi=300)
plt.close()

print(f"Generated comparison plots in {out_folder}")
