# Sampling transitions

## ... by direct simulation
These functions generate noise-induced transitions between an initial and final state.

```@docs
transition(sys::CoupledSDEs, x_i::State, x_f::State; kwargs...)
transitions(sys::CoupledSDEs, x_i::State, x_f::State, N=1; kwargs...)
```

## ... in pathspace
<!-- 
```@docs
langevinmcmc(sys::CoupledSDE, init; kwargs...)
stochastic_bridge(sys::CoupledSDE, Tphys::Float64, Δz::Float64)
symbolise_spde(sys::CoupledSDE)
langevinmcmc_spde(u, p, t)
``` -->