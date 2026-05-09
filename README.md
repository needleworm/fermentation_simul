# An ODE-Network-Based Chemical Complex System Simulator for Multi-Environment Fermentation

This repository contains the official implementation of the paper:

> **An ODE-Network-Based Chemical Complex System Simulator for Multi-Environment Fermentation**

The simulator models fermentation as a **chemical complex system** in which substrates, enzymes, yeast biomass, oxygen, pH, temperature, and metabolic products co-evolve through an interpretable ODE-network structure.

The framework supports multiple fermentation environments using a unified model skeleton:

- Bread dough fermentation
- Rice wine fermentation
- Beer wort fermentation
- Malt mash / wort fermentation

---

# Features

- ODE-network-based fermentation simulator
- Multi-environment fermentation presets
- Coupled saccharification–fermentation dynamics
- Interpretable state-transition structure
- Automated condition sweep generation
- Failure-mode discovery
- CSV-based experiment export
- ODE-network visualization tools

---

# Repository Structure

```text
.
├── fermentation_benchmark_full.py
├── visualize_ode_network.py
├── fermentation_benchmark_results/
├── ode_network_visualization/
└── README.md
```

---

# ODE-Network Structure

The simulator does not store the ODE network as an explicit graph object.

Instead:

- `STATE_NAMES` defines dynamic state variables
- `rhs()` defines flux equations and interactions
- `euler_step()` performs synchronous ODE updates

The network includes:

- substrate conversion
- enzymatic saccharification
- yeast uptake
- oxygen-dependent growth
- ethanol inhibition
- acidification
- CO₂ generation
- gas retention

The structural graph can be reconstructed using:

```bash
python visualize_ode_network.py
```

Generated outputs:

```text
ode_network_visualization/
├── ode_network_core.png
├── ode_network_full.png
├── ode_network_full_edges.csv
└── ode_network_full_nodes.csv
```

---

# Running the Simulator

## Generate Full Benchmark Sweep

```bash
python fermentation_benchmark_full.py
```

This automatically generates multiple fermentation conditions and exports CSV results.

Output directory:

```text
fermentation_benchmark_results/
```

Generated files:

```text
experiment_conditions.csv
timeseries.csv
summary.csv
environment_rankings.csv
column_dictionary.csv
manifest.csv
```

---

# Fermentation Environments

## Bread Dough

Focus:

- CO₂ productivity
- retained gas
- dough-like fermentation dynamics

## Rice Wine

Focus:

- simultaneous saccharification–fermentation
- ethanol yield
- pH dynamics

## Beer Wort

Focus:

- attenuation
- ethanol production
- residual sugar dynamics

## Malt Mash / Wort

Focus:

- enzymatic starch conversion
- mash temperature effects
- coupled saccharification–fermentation

---

# Validation

The simulator produces literature-aligned semi-quantitative trends, including:

- Rice-wine fermentation near 25°C producing approximately 11% v/v ethanol
- Malt-centered mash conditions near 62°C
- Bread-dough CO₂ productivity near 30°C
- Ale-like beer fermentation behavior

The simulator is intended as a:

- lightweight research scaffold
- interpretable ODE-network benchmark
- chemical complex system simulator

rather than a fully calibrated industrial production model.

---

# Citation

If you use this repository in your research, please cite:

```text
Byunghyun Ban, 2026, 
"An ODE-Network-Based Chemical Complex System Simulator for Multi-Environment Fermentation", Researchgate Preprint.

DOI:
https://doi.org/10.13140/RG.2.2.31603.28968
```

---

# Related Works

This repository extends previous studies on ODE-network-based chemical complex systems and agricultural nutrient interaction models:

- Machine learning approach to remove ion interference effect in agricultural nutrient solutions
- Nutrient Solution Management System for Smart Farms and Plant Factory
- ODE Network Model for Nonlinear and Complex Agricultural Nutrient Solution System
- Mathematical Model and Simulation for Nutrient-Plant Interaction Analysis
- The phenotype control kernel of a biomolecular regulatory network

---

# License

This repository is released for academic and research purposes.
