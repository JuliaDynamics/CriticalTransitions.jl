# Sampling transitions

## ... by direct simulation
These functions generate noise-induced transitions between an initial and final state.

```@docs
transition(sys::StochSystem, x_i::State, x_f::State; kwargs...)
transitions(sys::StochSystem, x_i::State, x_f::State, N=1; kwargs...)
```

## ... in pathspace
```@docs
langevinmcmc(sys::StochSystem, init; kwargs...)
```