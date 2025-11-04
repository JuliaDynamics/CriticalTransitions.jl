# CriticalTransitions.jl

*A Julia package for numerical investigation of noise- and rate-induced transitions in dynamical systems.*

This package provides a toolbox for dynamical systems under time-dependent forcing, with a focus on tipping phenomena and metastability. 
Building on [DynamicalSystems.jl](https://juliadynamics.github.io/DynamicalSystems.jl/stable/) and [DifferentialEquations.jl](https://diffeq.sciml.ai/stable/), the code adds functionality specifically to study stochastic and non-autonomous dynamics -- while keeping a familiar user interface and taking advantage of the powerful existing solvers in Julia.

## Purpose
The DynamicalSystems.jl ecosystem offers various functions to simulate and analyze dynamical systems in continuous time, centered around the `CoupledODEs` type that interfaces with the solvers of DifferentialEquations.jl. However, dedicated functionality for stochastic and/or explicitly time-dependent dynamics has largely been missing in this framework. [CriticalTransitions.jl](https://github.com/JuliaDynamics/CriticalTransitions.jl) seeks to fill this gap.

## Features
!!! tip "With this package, you can"
    - easily construct stochastic and nonautonomous dynamical systems
    - efficiently sample transition path ensembles
    - calculate minimum action paths and critical forcing rates
    - use a growing toolbox of tested and documented functions implementing concepts of large deviation theory, transition path theory, and rate-induced tipping

    ... and more, check out the [Examples](@ref) and [Manual](@ref) sections of these docs!

!!! ukw "Planned features"
    * **Rare event simulation**: importance sampling, AMS
    * **Quasipotentials**: Ordered line integral method (OLIM)
    * **Symbolic** differentiation of action functionals
    
    Missing a feature? [Open an issue](https://github.com/JuliaDynamics/CriticalTransitions.jl/issues) and tell us about it!

## People
Main developers:
- Reyk Börner (@reykboerner)
- Orjan Ameye (@oameye)
- Ryan Deeley (@ryandeeley)
- Raphael Römer (@raphael-roemer)

Thanks to George Datseris, Jeroen Wouters, Calvin Nesbitt, Tobias Grafke & Oliver Mehling.

This package got started in the EU-funded [CriticalEarth](https://www.criticalearth.eu) project.