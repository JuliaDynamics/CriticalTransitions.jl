# Simulating the system

We provide two main functions to simulate a [`StochSystem`](@ref) forward in time:
* [`relax`](@ref) for the system's deterministic evolution (in the absence of noise) based on the [`ODEProblem`](https://diffeq.sciml.ai/stable/types/ode_types/#SciMLBase.ODEProblem) of `DifferentialEquations.jl`
* [`simulate`](@ref) for stochastic dynamics based on the [`SDEProblem`](https://diffeq.sciml.ai/stable/types/sde_types/#SciMLBase.SDEProblem) of `DifferentialEquations.jl`

## Deterministic dynamics
```@docs
relax(sys::StochSystem, init::State; kwargs...)
```

## Stochastic dynamics

```@docs
simulate(sys::StochSystem, init::State; kwargs...)
```