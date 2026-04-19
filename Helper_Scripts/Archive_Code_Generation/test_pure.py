import numpy as np
from itertools import product
from qutip import *

def generate_combinations(K, N):
    all_combinations = product(range(N + 1), repeat=K)
    return [comb for comb in all_combinations if sum(comb) <= N]

def initialize_HB_state(N, d, D=None):
    D = N + 1 if D is None else D
    modes = d + 1
    n_list = [N // modes] * modes
    for i in range(N % modes):
        n_list[i] += 1
        
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
            b_k += (1.0 / np.sqrt(modes)) * (omega ** (j * k)) * a_ops[j]
        b_ops.append(b_k)
        
    state = tensor([basis(D, 0)] * modes)
    for k in range(modes):
        for _ in range(n_list[k]): state = b_ops[k] * state

    state = state.unit() 
    combinations = generate_combinations(d, N)

    coeffs = []
    for comb in combinations:
        n_0 = N - sum(comb)
        target_ket = tensor([basis(D, n_0)] + [basis(D, n_i) for n_i in comb])
        coeffs.append(target_ket.overlap(state))
        
    return np.array(coeffs), combinations

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

N = 3
d = 2
coeffs, combs = initialize_HB_state(N, d)
qfim_pure = calculate_QFIM_pure(coeffs, combs, d)
qfi_max = np.trace(qfim_pure)
crb = np.trace(np.linalg.pinv(qfim_pure))
print(f"Pure state calculation N={N}, d={d} -> QFI={qfi_max:.3f}, CRB={crb:.3f}")

# Compare to known values from notebook output for eta=1.0, p=1.0, N=3, d=2: QFI=10.667, CRB=0.500
