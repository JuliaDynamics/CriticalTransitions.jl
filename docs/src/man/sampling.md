# Sampling transitions

## ... by direct simulation
These functions generate noise-induced transitions between an initial and final state.

```@docs
transition(sys::StochSystem, x_i::State, x_f::State; kwargs...)
transitions(sys::StochSystem, x_i::State, x_f::State, N=1; kwargs...)
```

## ... in pathspace

```@docs
langevinmcmc_spde(u, p, t)
```

```@docs
symbolise_spde(sys::StochSystem)
```

```@docs
stochtolangevinmcmcstoch(sys::StochSystem, Tphys::Float64, Δz::Float64)
```

```@docs
langevinmcmc(sys::StochSystem, init; kwargs...)
```