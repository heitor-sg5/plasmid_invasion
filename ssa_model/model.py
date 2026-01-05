import numpy as np

y0 = [0.9, 0.1]

def analyze_ssa_results(results):
    for s, data in results.items():
        F_matrix = data['F']
        P_matrix = data['P']
        mean_F = F_matrix.mean(axis=0)
        std_F = F_matrix.std(axis=0)
        mean_P = P_matrix.mean(axis=0)
        std_P = P_matrix.std(axis=0)
        extinct_fraction = np.sum(P_matrix[:, -1] == 0) / P_matrix.shape[0]
        print(f"s = {s:.2f}")
        print(f"Mean F (final) = {mean_F[-1]:.4f} ± {std_F[-1]:.4f}")
        print(f"Mean P (final) = {mean_P[-1]:.4f} ± {std_P[-1]:.4f}")
        print(f"Fraction of runs where plasmids go extinct: {extinct_fraction:.2f}\n")

def gillespie_ssa(p, TMAX, y0):
    F = int(y0[0] * p['K'])
    P = int(y0[1] * p['K'])
    t = 0.0

    K = p['K']
    r = p['r']
    mu = p['mu']
    delta = p['delta']
    beta = p['beta']
    c = p['c']
    s = p['s']

    rF = r * (1 - s)
    rP = r * (1 - c)
    K_div = 1.0 / K
    beta_div = beta * K_div

    times = [t]
    Fs = [F]
    Ps = [P]

    while t < TMAX and (F + P) > 0:
        N = F + P
        NFactor = 1 - N * K_div

        wF = max(0.0, rF * NFactor)
        wP = max(0.0, rP * NFactor)
        a0 = wF * F + wP * P + mu * (F + P) + delta * P + beta_div * F * P

        if a0 <= 0:
            break

        tau = np.random.exponential(1 / a0)
        t += tau

        r_val = np.random.rand() * a0

        c0 = wF * F
        if r_val < c0:
            F += 1
        else:
            c1 = c0 + wP * P
            if r_val < c1:
                P += 1
            else:
                c2 = c1 + mu * F
                if r_val < c2:
                    F -= 1
                else:
                    c3 = c2 + mu * P
                    if r_val < c3:
                        P -= 1
                    else:
                        c4 = c3 + delta * P
                        if r_val < c4:
                            P -= 1
                            F += 1
                        else:
                            F -= 1
                            P += 1
        F = max(F, 0)
        P = max(P, 0)

        times.append(t)
        Fs.append(F)
        Ps.append(P)

    return np.array(times), np.array(Fs), np.array(Ps)

def run_multiple_ssa(base_params, TMAX, runs):
    results = {}
    s_values=[0.0, 0.2, 0.6]
    t_grid = np.linspace(0, TMAX, 500)
    base_params['K'] = base_params['K'] * 1000

    for s in s_values:
        Fs_runs = []
        Ps_runs = []

        for _ in range(runs):
            p = base_params.copy()
            p['s'] = s
            t, F, P = gillespie_ssa(p, TMAX, y0)

            F_interp = np.interp(t_grid, t, F)
            P_interp = np.interp(t_grid, t, P)

            Fs_runs.append(F_interp / p['K'])
            Ps_runs.append(P_interp / p['K'])

        results[s] = {
            't': t_grid,
            'F': np.array(Fs_runs),
            'P': np.array(Ps_runs)
        }

    return results