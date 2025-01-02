# Large deviation theory


## Mimumum action path
```@docs
CriticalTransitions.MinimumActionPath
```

## String method
The string method is a technique for finding transition paths between two states in a dynamical system. The method represents the path as a "string" of points that connects the states and evolves it to minimize the drift along the path. The resulating  tangent path is parrallel to the drift of the system, i.e., the string method computes the heteroclinic orbit. For non-gradient systems (detailed -balance is broken), the heteroclinic orbit differs from the transition path, it does correctly predict, it correctly captures the deterministic dynamics from the saddle point onward ("downhill" portion of the path).
```@docs
string_method
```

## Minimum action methods
The minimum action method (MAM) is a technique for calculating the most probable transition path between two (meta)stable states in a stochastic dynamical system. In the limit of small noise, this path corresponds to the minimizer of an action functional. The action functional typically takes into account both the deterministic drift and the noise intensity of the system. By discretizing this path and using optimization techniques, MAM finds the trajectory that requires the least "effort" to transition between states in phase space.

### Minimum action method (MAM)
Minimization of the action functioal using the optimization algorithm of `Optimization.jl`.

```@docs
min_action_method
```

### Geometric minimum action method (gMAM)
Minimization of the geometric action following
[Heymann and Vanden-Eijnden, PRL (2008)](https://link.aps.org/doi/10.1103/PhysRevLett.100.140601).
The gMAM reformulates MAM to avoid double optimisation of both the action and the transition time. It achieves this by using a geometric action functional that is independent of the time parametrization of the path. This reparametrization invariance makes the method more robust and computationally efficient, particularly for systems with metastable states separated by large barriers.
```@docs
geometric_min_action_method
```

### Simple Geometric minimum action method (sgMAM)
Simplified minimization of the geometric action following
[Grafke et al. (2017)](https://doi.org/10.1007/978-1-4939-6969-2_2).
The simple gMAM reduces the complexity of the original gMAM by requiring only first-order derivatives of the underlying Hamiltonian optimisation formulation. This simplifies the numerical treatment and the computational complexity.

The implementation below perform a constrained gradient descent where it assumes an autonomous system with additive Gaussian noise.
```@docs
sgmam
SgmamSystem
```

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