import os
import matplotlib.pyplot as plt

def ensure_figures_dir():
    fig_dir = os.path.join(os.path.dirname(__file__), 'figures')
    if not os.path.exists(fig_dir):
        os.makedirs(fig_dir)
    return fig_dir

def plot_time_series(results):
    fig_dir = ensure_figures_dir()
    s_values = sorted(results.keys())

    fig, axs = plt.subplots(1, len(s_values), figsize=(18,5))
    for i, s in enumerate(s_values):
        sol = results[s]
        F, P = sol.y
        axs[i].plot(sol.t, F, label='F')
        axs[i].plot(sol.t, P, label='P')
        axs[i].set_title(f'Time Series at s={s}')
        axs[i].set_xlabel('Time')
        axs[i].set_ylabel('Population size')
        axs[i].legend()
    plt.tight_layout()
    plt.show()
    fig_path = os.path.join(fig_dir, 'time_series.png')
    fig.savefig(fig_path)
    print(f"Time series figure saved to {fig_path}")

def plot_beta_heatmap(heatmap, beta_values, cost_values, beta_sweep, pf_sweep):
    fig_dir = ensure_figures_dir()
    fig, axs = plt.subplots(1,2, figsize=(12,5))
    im = axs[0].imshow(heatmap, aspect='auto', origin='lower', 
                       extent=[beta_values[0], beta_values[-1],
                       cost_values[0], cost_values[-1]], cmap='viridis')
    axs[0].set_xlabel('Transfer rate, β')
    axs[0].set_ylabel('Plasmid cost, c')
    axs[0].set_title('Plasmid Fraction Heatmap')
    fig.colorbar(im, ax=axs[0])

    axs[1].plot(beta_sweep, pf_sweep, marker='o')
    axs[1].set_xlabel('Transfer rate, β')
    axs[1].set_ylabel('Plasmid frequency')
    axs[1].set_title('β Bifurcation Plot')
    axs[1].grid(True)

    plt.tight_layout()
    plt.show()
    fig_path = os.path.join(fig_dir, 'beta_heatmap.png')
    fig.savefig(fig_path)
    print(f"Beta heatmap figure saved to {fig_path}")
