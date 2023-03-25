# CriticalTransitions.jl

*A Julia package for the numerical investigation of stochastic dynamical systems.*

Building on [`DifferentialEquations.jl`](https://diffeq.sciml.ai/stable/) and [`DynamicalSystems.jl`](https://juliadynamics.github.io/DynamicalSystems.jl/stable/), this package aims to provide a toolbox for stochastic dynamical systems, with a focus on tipping phenomena and metastability.

!!! compat "Compatibility"
    This is the latest working version (`v0.1.1`) compatible with `DynamicalSystems.jl v2` (specifically `v2.3.2`). Our next release will contain breaking changes in order to introduce a modified type structure and to be compatible with `DynamicalSystems.jl v3`.


!!! info "Current features"
    * **Stability analysis**: Fixed points, linear stability, basins of attraction, edge tracking
    * **Stochastic simulation**: Gaussian noise, uncorrelated and correlated, additive and multiplicative
    * **Transition path sampling**: Ensemble sampling by direct simulation and Pathspace Langevin MCMC
    * **Large deviation theory**: Action functionals and minimization algorithms (MAM, gMAM)

!!! ukw "Planned features"
    * **Rare event simulation**: importance sampling, AMS
    * **Quasipotentials**: Ordered line integral method (OLIM)
    * **Rate-induced tipping** tools
    * ...?


Developers: Reyk Börner, Ryan Deeley and Raphael Römer

Further contributors: Jeroen Wouters, Oliver Mehling

This work is part of the [CriticalEarth](https://criticalearth.eu) project.