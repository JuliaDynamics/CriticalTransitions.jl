# CriticalTransitions.jl

*A Julia package for the numerical investigation of critical transitions in stochastic dynamical systems.*

Building on [`DifferentialEquations.jl`](https://diffeq.sciml.ai/stable/) and [`DynamicalSystems.jl`](https://juliadynamics.github.io/DynamicalSystems.jl/stable/), this newly developing package aims to provide a toolbox for stochastic dynamical systems, with a focus on tipping phenomena and metastability.

!!! info "Current features"
    * **Stability analysis**: fixed points, linear stability, basins of attraction, basin boundary
    * **Stochastic simulation**: Gaussian noise, uncorrelated and correlated, additive and multiplicative
    * **Transition path sampling**: parallelized ensemble sampling
    * **Large deviation theory**: Freidlin-Wentzell and Onsager-Machlup action calculation

!!! ukw "Planned features"
    * **Rare event simulation**: importance sampling, AMS
    * **Large deviations**: instantons and quasipotentials
    * **Edge tracking**
    * **Langevin MCMC path sampling**
    * ...?


Developers: Reyk Börner, Ryan Deeley, and Raphael Römer

This work is part of the [CriticalEarth](https://criticalearth.eu) project.