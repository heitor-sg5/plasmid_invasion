import time
from ode_model import model as ode_model, charts as ode_charts
from ssa_model import model as ssa_model, charts as ssa_charts

def main():
    print("Select simulation type:")
    print("1 - ODE model")
    print("2 - SSA model")
    choice = input("Enter 1 or 2: ").strip()

    if choice == '1':
        start_time = time.time()
        p = ode_model.base_params.copy()
        beta_crit = ode_model.get_beta_crit(p)
        s_crit = ode_model.get_s_crit(p)

        print("\nAnalytical thresholds (base parameters):")
        print(f"Critical β = {beta_crit:.4f}")
        print(f"Critical s = {s_crit:.4f}\n")

        ts_results = ode_model.run_time_series()
        print(f"Time series simulation completed in {time.time() - start_time:.2f} seconds.\n")

        start_time = time.time()
        ode_charts.plot_time_series(ts_results)
        heatmap, beta_values, cost_values, beta_sweep, pf_sweep = ode_model.run_beta_sweep()
        print(f"Beta sweep simulation completed in {time.time() - start_time:.2f} seconds.\n")
        
        ode_charts.plot_beta_heatmap(heatmap, beta_values, cost_values, beta_sweep, pf_sweep)
    elif choice == '2':
        start_time = time.time()
        results = ssa_model.run_multiple_ssa()
        print(f"SSA simulation completed in {time.time() - start_time:.2f} seconds.\n")
        ssa_charts.plot_ssa_trajectories(results)
        
    else:
        print("Invalid choice.")

if __name__ == "__main__":
    main()