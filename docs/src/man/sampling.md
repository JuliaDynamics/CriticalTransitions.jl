# Sampling transitions

## ... by direct simulation
These functions generate noise-induced transitions between an initial and final state.

```@docs
transition(sys::CoupledSDEs, x_i, x_f; kwargs...)
transitions(sys::CoupledSDEs, x_i, x_f, N=1; kwargs...)
```

## ... in pathspace