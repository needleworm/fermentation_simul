ODE-network visualization outputs
=================================

The ODE-network is implicit in fermentation_benchmark_full.py: STATE_NAMES defines dynamic nodes, rate_factors()/stress_index() define regulatory factors, rhs() defines fluxes, and euler_step() performs synchronous updates.

Recommended manuscript figure: ode_network_core.svg or ode_network_core.png.
Recommended supplementary figure: ode_network_full.svg or ode_network_full.png.
Edge-list CSVs are included for transparent reporting.
