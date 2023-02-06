# Large deviation theory

## Action functionals

### Freidlin-Wentzell action
```@docs
fw_action(sys::StochSystem, path, time; kwargs...)
```

### Geometric Freidlin-Wentzell action
```@docs
geometric_action(sys::StochSystem, path, arclength=1; kwargs...)
```

### Onsager-Machlup action
```@docs
om_action(sys::StochSystem, path, time; kwargs...)
```

## Minimum action paths
We provide the following two methods to calculate *instantons*, or minimum action paths,
between two states of a `StochSystem`.

### Minimum action method (MAM)
Minimization of the Freidlin-Wentzell action using the L-BFGS algorithm of `Optim.jl`.

```@docs
mam(sys::StochSystem, x_i::State, x_f::State, N::Int, T::Float64; kwargs...)
```

### Geometric minimum action method (gMAM)
Minimization of the geometric action following
[Heymann and Vanden-Eijnden, PRL (2008)](https://link.aps.org/doi/10.1103/PhysRevLett.100.140601).

```@docs
gmam(sys::StochSystem, x_i::State, x_f::State, arclength=1; kwargs...)
```