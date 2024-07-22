# Large deviation theory

## Action functionals

### Freidlin-Wentzell action
```@docs
fw_action
```

### Geometric Freidlin-Wentzell action
```@docs
geometric_action
```

### Onsager-Machlup action
```@docs
om_action
```

For convenience, a general [`action`](@ref) function is available where the type of functional is set as an argument:

```@docs
action
```

## Minimum action paths
We provide the following two methods to calculate *instantons*, or minimum action paths,
between two states of a `CoupledSDEs`.

### Minimum action method (MAM)
Minimization of the Freidlin-Wentzell action using the L-BFGS algorithm of `Optim.jl`.

```@docs
min_action_method
```

### Geometric minimum action method (gMAM)
Minimization of the geometric action following
[Heymann and Vanden-Eijnden, PRL (2008)](https://link.aps.org/doi/10.1103/PhysRevLett.100.140601).

```@docs
geometric_min_action_method
```