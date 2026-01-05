import os
import matplotlib.pyplot as plt

def ensure_figures_dir():
    fig_dir = os.path.join(os.path.dirname(__file__), 'figures')
    if not os.path.exists(fig_dir):
        os.makedirs(fig_dir)
    return fig_dir

def plot_ssa_trajectories(results):
    fig_dir = ensure_figures_dir()
    s_values = sorted(results.keys())
    fig, axs = plt.subplots(1, 3, figsize=(18, 5), sharey=True)

    for i, s in enumerate(s_values):
        data = results[s]
        t = data['t']
        F = data['F']
        P = data['P']

        for k in range(F.shape[0]):
            axs[i].plot(t, F[k], color='tab:blue', alpha=0.2)
            axs[i].plot(t, P[k], color='tab:orange', alpha=0.2)
        axs[i].plot(t, F.mean(axis=0), color='tab:blue', linewidth=3, label='F (mean)')
        axs[i].plot(t, P.mean(axis=0), color='tab:orange', linewidth=3, label='P (mean)')
        axs[i].set_title(f'Time Series at s={s}')
        axs[i].set_xlabel('Time')
        if i == 0:
            axs[i].set_ylabel('Population size')
        axs[i].legend()

    plt.tight_layout()
    plt.show()
    fig_path = os.path.join(fig_dir, 'ssa_trajectories.png')
    fig.savefig(fig_path)
    print(f"SSA trajectories saved to {fig_path}")