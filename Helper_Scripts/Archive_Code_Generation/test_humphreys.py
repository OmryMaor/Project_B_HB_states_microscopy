import numpy as np
from qutip import *
from itertools import product

def initialize_true_HB_state(n, d, D=None):
    # Total photons N = n * (d+1)
    N = n * (d+1)
    D = N + 1 if D is None else D
    modes = d + 1
    
    # Input state is |n, n, n, ... n> in modes
    omega = np.exp(2j * np.pi / modes)
    a_ops = []
    for i in range(modes):
        ops = [qeye(D)] * modes
        ops[i] = create(D)
        a_ops.append(tensor(ops))
        
    b_ops = []
    for k in range(modes):
        b_k = 0
        for j in range(modes):
            # QFT matrix
            b_k += (1.0 / np.sqrt(modes)) * (omega ** (j * k)) * a_ops[j]
        b_ops.append(b_k)
        
    state = tensor([basis(D, 0)] * modes)
    for k in range(modes):
        for _ in range(n): 
            state = b_ops[k] * state
    state = state.unit() 
    
    # Generate combinations
    all_combinations = product(range(N + 1), repeat=d)
    combinations = [comb for comb in all_combinations if sum(comb) <= N]

    coeffs = []
    for comb in combinations:
        n_0 = N - sum(comb)
        target_ket = tensor([basis(D, n_0)] + [basis(D, n_i) for n_i in comb])
        coeffs.append(target_ket.overlap(state))
        
    return np.array(coeffs), combinations, N

def calculate_QFIM_pure(coeffs, combs, d):
    K = d
    qfim = np.zeros((K, K))
    n_exp = np.zeros(K)
    for idx, c in enumerate(coeffs):
        p = np.abs(c)**2
        for a in range(K):
            n_exp[a] += p * combs[idx][a]
            
    n_corr = np.zeros((K, K))
    for idx, c in enumerate(coeffs):
        p = np.abs(c)**2
        for a in range(K):
            for b in range(K):
                n_corr[a, b] += p * combs[idx][a] * combs[idx][b]
                
    for a in range(K):
        for b in range(K):
            qfim[a, b] = 4 * (n_corr[a, b] - n_exp[a] * n_exp[b])
    return qfim

print("Testing d=3 with n=1 (Total N=4)")
c, combs, N_tot = initialize_true_HB_state(1, 3)
qf = calculate_QFIM_pure(c, combs, 3)
crb_hb = np.trace(np.linalg.pinv(qf))
print(f"HB(1, 3) CRB = {crb_hb:.3f}")

print("Testing d=3 with n=2 (Total N=8)")
c, combs, N_tot = initialize_true_HB_state(2, 3)
qf = calculate_QFIM_pure(c, combs, 3)
crb_hb = np.trace(np.linalg.pinv(qf))
print(f"HB(2, 3) CRB = {crb_hb:.3f}")

print("Testing d=4 with n=1 (Total N=5)")
c, combs, N_tot = initialize_true_HB_state(1, 4)
qf = calculate_QFIM_pure(c, combs, 4)
crb_hb = np.trace(np.linalg.pinv(qf))
print(f"HB(1, 4) CRB = {crb_hb:.3f}")
