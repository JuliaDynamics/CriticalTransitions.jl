# Simulating the system

## SDE integration

```@docs
simulate(sys::StochSystem, init::State; kwargs...)
relax(sys::StochSystem, init::State; kwargs...)
```

## Noise-induced transitions

```@docs
transition(sys::StochSystem, x_i::State, x_f::State; kwargs...)
transitions(sys::StochSystem, x_i::State, x_f::State, N=1; kwargs...)
```
