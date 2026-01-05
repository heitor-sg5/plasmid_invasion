# Plasmid Dynamics Models

## Overview

Plasmid invasion dynamics are studied using two complementary modelling approaches: 
- Deterministic ordinary differential equations (ODEs) 
- Stochastic simulation algorithm (SSA) 

Both models track the population dynamics of plasmid-free bacteria (F) and plasmid-bearing bacteria (P), incorporating: 
- Density-dependent population growth
- Plasmid transfer via conjugation 
- Segregational plasmid loss
- Plasmid-associated fitness costs 
- Selective pressure 

The deterministic model captures mean-field population behaviour, while the stochastic model resolves demographic noise and rare-event dynamics that arise when populations are small. Together, these approaches allow comparison between average behaviour and stochastic variability in plasmid persistence and invasion.

---

## Project Structure

```
├── main.py
├── parameters/
│   ├── params.json
├── ode_model/
│   ├── model.py
│   ├── charts.py
│   └── figures/
├── ssa_model/
│   ├── model.py
│   ├── charts.py
│   └── figures/
├── docs/
│   ├── appendix_A.pdf
│   ├── appendix_B.pdf
│   ├── appendix_C.pdf
│   └── analysis.pdf
├── requirements.txt
└── README.md
```

---

## Installation

Prerequisites:
- Python 3.8 or higher
- NumPy
- SciPy
- Matplotlib

1. Close the repository:

```bash
git clone https://github.com/heitor-sg5/plasmid_invasion.git
cd plasmid-dynamics-models
```

2. Install Python dependencies:

```bash
pip install -r requirements.txt
```

---

## Usage

Run the program from the project root:

```bash
python main.py
```

This will:

1. load the parameters from `parameters/params.json`
2. Run deterministic ODE time-series simulations
3. Perform a β-cost parameter sweep
4. Compute analytical threshold conditions
5. Generate and save figures to `ode_model/figures`

For the stochastic SSA simulation:

```bash
python main.py --ssa
```

Arguments may be customized, for instance:

```bash
python main.py --ssa --runs 100 --TMAX 1000
```

Available CL options:

- `--params`: path to parameter JSON file
- `--TMAX`: maximum simulation time
- `--ode`: run ODE model
- `--ssa`: run SSA model
- `--runs`: number of SSA simulations

---

## Documentation

A detailed analysis of the model, including derivations of critical conjugation thresholds and interpretation of simulation results, are provided in the PDF files found in `docs/`.