# Simulating the system

We provide two main functions to simulate a [`StochSystem`](@ref) forward in time:
* [`relax`] for the system's deterministic evolution (in the absence of noise) based on the [`ODEProblem`] of `DifferentialEquations.jl`
* [`simulate`] for stochastic dynamics based on the [`SDEProblem`] of `DifferentialEquations.jl`

## Deterministic dynamics
```@docs
relax(sys::StochSystem, init::State; kwargs...)
```

## Stochastic dynamics

```@docs
simulate(sys::StochSystem, init::State; kwargs...)
```