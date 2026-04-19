import numpy as np
import pandas as pd
from qutip import *
from itertools import product
from scipy.optimize import minimize
import os
import csv

def generate_combinations(K, N):
    all_combinations = product(range(N + 1), repeat=K)
    return [comb for comb in all_combinations if sum(comb) <= N]

def diagonalize(rho):
    return rho.eigenstates()

def compute_rho_derivative(rho, vals, vecs, generators, a):
    n_a = generators[a]
    U = Qobj(np.hstack([v.full() for v in vecs]), dims=[rho.dims[0], rho.dims[0]])
    n_a_eigen = (U.dag() * n_a * U).full()
    diffs = vals[np.newaxis, :] - vals[:, np.newaxis]
    return -1j * diffs * n_a_eigen

def calculate_QFIM(rho, vals, vecs, generators):
    dim = len(vals)
    K = len(generators)
    qfim = np.zeros((K, K))
    derivatives_eigen = [compute_rho_derivative(rho, vals, vecs, generators, a) for a in range(K)]

    for a in range(K):
        for b in range(a, K):
            term_sum = 0
            for n in range(dim):
                for m in range(dim):
                    v_sum = vals[n] + vals[m]
                    if v_sum > 1e-14:
                        elem_a = derivatives_eigen[a][n, m]
                        elem_b = derivatives_eigen[b][m, n]
                        term_sum += (2.0 / v_sum) * np.real(elem_a * elem_b)
            qfim[a, b] = term_sum
            qfim[b, a] = term_sum
    return qfim

# ======== PRE-CACHED MATRICES FOR SPEED ========
def get_cached_V_free(eta_vec, N, K):
    D = N + 1
    V_total = tensor([qeye(D)] * (2 * K))
    for i in range(K):
        theta_i = np.arccos(np.sqrt(eta_vec[i]))
        ops_s = [qeye(D)] * (2 * K); ops_s[i] = destroy(D)
        ops_e = [qeye(D)] * (2 * K); ops_e[K + i] = destroy(D)
        h_int = tensor(ops_s).dag() * tensor(ops_e) - tensor(ops_s) * tensor(ops_e).dag()
        V_total = (theta_i * h_int).expm() * V_total
    return V_total

def get_cached_V_int(eta_vec, N, K):
    D = N + 1
    a_e = tensor([qeye(D)] * K + [destroy(D)])
    h_tot = tensor([qeye(D)] * (K + 1)) * 0.0
    for i in range(K):
        theta_i = np.arccos(np.sqrt(eta_vec[i]))
        ops_s = [qeye(D)] * (K + 1); ops_s[i] = destroy(D)
        h_tot += theta_i * (tensor(ops_s).dag() * a_e - tensor(ops_s) * a_e.dag())
    return h_tot.expm()

def apply_N_free(initial_signal_state, V_total, N, K):
    D = N + 1
    rho_joint = ket2dm(tensor(initial_signal_state, tensor([basis(D, 0)] * K)))
    return (V_total * rho_joint * V_total.dag()).ptrace(list(range(K)))

def apply_N_int(initial_signal_state, V_total, N, K):
    D = N + 1
    rho_joint = ket2dm(tensor(initial_signal_state, basis(D, 0)))
    return (V_total * rho_joint * V_total.dag()).ptrace(list(range(K)))

def apply_N_tot(initial_ket, V_free, V_int, p, N, K):
    if p == 1.0: return apply_N_free(initial_ket, V_free, N, K)
    if p == 0.0: return apply_N_int(initial_ket, V_int, N, K)
    return p * apply_N_free(initial_ket, V_free, N, K) + (1.0 - p) * apply_N_int(initial_ket, V_int, N, K)

def optimize_initial_state_CRB(N, K, eta_vec, p):
    D = N + 1
    combinations = generate_combinations(K, N)
    num_coeffs = len(combinations)

    generators = []
    for i in range(K):
        ops = [qeye(D)] * K; ops[i] = num(D)
        generators.append(tensor(ops))

    basis_kets = [tensor([basis(D, n) for n in comb]) for comb in combinations]

    # Extreme Performance Gain: Cache physical mappings
    V_free = get_cached_V_free(eta_vec, N, K)
    V_int = get_cached_V_int(eta_vec, N, K)

    def objective(coeffs):
        norm = np.linalg.norm(coeffs)
        if norm < 1e-10: return 1000.0
        
        initial_ket = sum((coeffs[i] / norm) * basis_kets[i] for i in range(num_coeffs))
        rho_out = apply_N_tot(initial_ket, V_free, V_int, p, N, K)
        vals, vecs = diagonalize(rho_out)
        
        qfim = calculate_QFIM(rho_out, vals, vecs, generators)
        
        if np.abs(np.linalg.det(qfim)) < 1e-10:
            return 1000.0
            
        return np.trace(np.linalg.inv(qfim))

    start_coeffs = np.array([1.0 / np.sqrt(num_coeffs)] * num_coeffs)
    cons = ({'type': 'eq', 'fun': lambda x: np.sum(x**2) - 1})
    bnds = [(0, 1) for _ in range(num_coeffs)]
    
    # We strictly use a single start due to high resolution and cache constraints
    res = minimize(objective, start_coeffs, method='SLSQP', bounds=bnds, constraints=cons, options={'maxiter': 50})

    opt_coeffs = res.x
    norm = np.linalg.norm(opt_coeffs)
    initial_ket = sum((opt_coeffs[i] / norm) * basis_kets[i] for i in range(num_coeffs))

    rho_out = apply_N_tot(initial_ket, V_free, V_int, p, N, K)
    vals, vecs = diagonalize(rho_out)
    qfim = calculate_QFIM(rho_out, vals, vecs, generators)

    max_qfi = np.trace(qfim)
    try:
        if np.abs(np.linalg.det(qfim)) < 1e-10:
             total_var = np.inf
        else:
             total_var = np.trace(np.linalg.inv(qfim))
    except:
        total_var = np.inf

    return res.x, combinations, max_qfi, total_var

# ================= EXECUTION =================
if __name__ == "__main__":
    N_param, K_param = 3, 2
    eta_grid = [0.8, 0.9, 0.95, 1.0]
    p_grid = [round(x, 1) for x in np.linspace(0.0, 1.0, 11)]
    
    out_dir = "Project_A_files"
    os.makedirs(out_dir, exist_ok=True)
    db_file = os.path.join(out_dir, "database_projectA_CRB.csv")
    
    # Delete file if exists to prevent messy partial appends
    if os.path.isfile(db_file): os.remove(db_file)
    
    with open(db_file, mode='a', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['N', 'K', 'eta_vec', 'p', 'max_qfi', 'total_variance', 'coeffs', 'combinations'])
            
        print("Starting Accelerated CRB Database Generation...")
        for eta in eta_grid:
            eta_vec = [eta, eta]
            for p in p_grid:
                print(f"Sweeping eta={eta}, p={p}...")
                coeffs, combs, qfi, var = optimize_initial_state_CRB(N_param, K_param, eta_vec, p)
                
                writer.writerow([N_param, K_param, str(eta_vec), p, qfi, var,
                                 ",".join(map(str, coeffs)),
                                 ";".join([",".join(map(str, c)) for c in combs])])
                f.flush()
    print("Database finalized!")
