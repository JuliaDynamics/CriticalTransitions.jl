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

For convenience, a general [`action`](@ref) function is available where the type of functional is set as an argument:

```@docs
action(sys::StochSystem, path::Matrix, time, functional; kwargs...)
```

## Minimum action paths
We provide the following two methods to calculate *instantons*, or minimum action paths,
between two states of a `StochSystem`.

### Minimum action method (MAM)
Minimization of the Freidlin-Wentzell action using the L-BFGS algorithm of `Optim.jl`.

```@docs
min_action_method(sys::StochSystem, x_i::State, x_f::State, N::Int, T::Float64; kwargs...)
```

### Geometric minimum action method (gMAM)
Minimization of the geometric action following
[Heymann and Vanden-Eijnden, PRL (2008)](https://link.aps.org/doi/10.1103/PhysRevLett.100.140601).

```@docs
geometric_min_action_method(sys::StochSystem, x_i::State, x_f::State, arclength=1; kwargs...)
```