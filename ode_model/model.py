import numpy as np
from scipy.integrate import solve_ivp

y0 = [0.99, 0.01]

def ode_model(t, y, p):
    F, P = y
    
    rF = p['r'] * (1 - p['s'])
    rP = p['r'] * (1 - p['c'])

    N = F + P
    wF = rF * (1 - N / p['K'])
    wP = rP * (1 - N / p['K'])

    mF = p['mu']
    mP = p['mu']

    dF = wF * F - mF * F + p['delta'] * P - p['beta'] * F * P
    dP = wP * P - mP * P - p['delta'] * P + p['beta'] * F * P

    return [dF, dP]

def get_beta_crit(p):
    num = p['delta'] + p['mu'] * (p['c'] / (1 - p['s']))
    denom = p['K'] * (1 - (p['mu'] / (p['r'] * (1 - p['s']))))
    return num / denom

def get_beta_cost(p):
    term1 = 1 / p['K']
    term2 = p['delta'] / (p['c'] * (1 - (p['mu'] / (p['r'] * (1 - p['s'])))))
    term3 = p['mu'] / ((1 - p['s'] * (1 - (p['mu'] / (p['r'] * (1 - p['s']))))))
    return term1 * (term2 + term3)

def get_beta_delta(p):
    term1 = 1 / (p['K'] * (1 - p['mu'] / (p['r'] * (1 - p['s']))))
    term2 = 1 + ((p['mu'] + p['c']) / (p['delta'] * (1 - p['s'])))
    return term1 * term2

def analyse_ode_results(base_params):
    p = base_params.copy()
    beta_crit = get_beta_crit(p)
    beta_cost = get_beta_cost(p)
    beta_delta = get_beta_delta(p)

    print("\nAnalytical thresholds (base parameters):")
    print(f"Critical β = {beta_crit:.4f}")
    print(f"Critical β / c = {beta_cost:.4f}")
    print(f"Critical β / δ = {beta_delta:.4f}")

def run_time_series(base_params, TMAX):
    results = {}
    s_values=[0.0, 0.2, 0.6]
    t_eval = np.linspace(0, TMAX, 500)
    for s in s_values:
        params = base_params.copy()
        params['s'] = s
        sol = solve_ivp(lambda t, y: ode_model(t, y, params), [0, TMAX], y0, t_eval=t_eval)
        results[s] = sol
    return results

def run_beta_sweep(base_params, TMAX):
    beta_values = np.linspace(0, 0.05, 20)
    cost_values = np.linspace(0, 0.2, 20)

    heatmap = np.zeros((len(cost_values), len(beta_values)))
    beta_sweep = []
    pf_sweep = []

    for i, c in enumerate(cost_values):
        for j, beta in enumerate(beta_values):
            params = base_params.copy()
            params['c'] = c
            params['beta'] = beta
            sol = solve_ivp(lambda t, y: ode_model(t, y, params), [0, TMAX], y0)
            F_final, P_final = sol.y[:, -1]
            heatmap[i, j] = P_final / (F_final + P_final)

    base_c = base_params['c']
    for beta in beta_values:
        params = base_params.copy()
        params['beta'] = beta
        params['c'] = base_c
        sol = solve_ivp(lambda t, y: ode_model(t, y, params), [0, TMAX], y0)
        F_final, P_final = sol.y[:, -1]
        beta_sweep.append(beta)
        pf_sweep.append(P_final / (F_final + P_final))

    return heatmap, beta_values, cost_values, beta_sweep, pf_sweep