# CriticalTransitions.jl

*A Julia package for numerical investigation of noise- and rate-induced transitions in dynamical systems.*

Building on [DynamicalSystems.jl](https://juliadynamics.github.io/DynamicalSystems.jl/stable/) and [DifferentialEquations.jl](https://diffeq.sciml.ai/stable/), this package provides a toolbox for dynamical systems under time-dependent forcing, with a focus on tipping phenomena and metastability.

![CT.jl infographic](./assets/CTjl_structure_v0.3_small.jpeg)

## Purpose
The `DynamicalSystems.jl` ecosystem offers various functions to simulate and analyze dynamical systems in continuous time, centered around the `CoupledODEs` type that interfaces with the powerful solvers of `DifferentialEquations.jl`. However, dedicated functionality for stochastic and/or explicitly time-dependent dynamics has largely been missing in this framework. `CriticalTransitions.jl` seeks to fill this gap: with this package, you can
- easily construct stochastic and nonautonomous dynamical systems
- efficiently generate trajectories, including paralellized ensemble simulation and transition path sampling
- use a growing toolbox of tested and documented functions implementing conecpts of large deviation theory, transition path theory, and rate-induced tipping

## Features

!!! info "Current features"
    * **Stochastic simulation** made easy: Gaussian noise, uncorrelated and correlated, additive and multiplicative
    * **Transition path sampling**: Parallelized ensemble rejection sampling
    * **Large deviation theory** tools: Action functionals and minimization algorithms (MAM, gMAM)

!!! ukw "Planned features"
    * **Rare event simulation**: importance sampling, AMS
    * **Quasipotentials**: Ordered line integral method (OLIM)
    * **Rate-induced tipping** tools
    * **Symbolic** differentiation of action functionals
    * ...?

## People
Main developers: Reyk Börner, Ryan Deeley, Raphael Römer and Orjan Ameye

Thanks to George Datseris, Jeroen Wouters, Calvin Nesbitt, Tobias Grafke, and Oliver Mehling.

This package got started in the EU-funded [CriticalEarth](https://www.criticalearth.eu) project.