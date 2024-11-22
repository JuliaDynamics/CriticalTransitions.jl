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
between two states of a `CoupledSDEs` system.

### Minimum action method (MAM)
Minimization of the Freidlin-Wentzell action using the L-BFGS algorithm of `Optim.jl`.

```@docs
min_action_method
```

### Geometric minimum action method (gMAM)
Minimization of the geometric action following
[Heymann and Vanden-Eijnden, PRL (2008)](https://link.aps.org/doi/10.1103/PhysRevLett.100.140601).
The gMAM reformulates MAM to avoid numerical stiffness by reparametrizing
the path in terms of arc length. This ensures an even distribution of points along the path,
enhancing numerical stability and accuracy. By solving the reparametrized integral,
gMAM accurately captures the geometry of the action functional and is well-suited for systems
with complex energy landscapes or intricate dynamics.

```@docs
geometric_min_action_method
```

### Simple Geometric minimum action method (sgMAM)
Simplified minimization of the geometric action following
[Heymann and Vanden-Eijnden, PRL (2008)](https://doi.org/10.1007/978-1-4939-6969-2_2).
The sgMAM is a streamlined version of the gMAM that simplifies the computation by avoiding
explicit reparametrization of the path. Instead, it introduces an implicit reparametrization
by rewriting the problem as a constrained optimization.

The implementation below perform a constrained gradient descent where it assumes an
autonomous system with additive Gaussian noise.
```@docs
SgmamSystem
sgmam
```