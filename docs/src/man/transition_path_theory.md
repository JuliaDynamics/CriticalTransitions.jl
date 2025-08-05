# Transition Path Theory

Transition Path Theory (TPT) provides a framework for analyzing rare transition events between metastable states in complex systems. This module implements TPT calculations for two-dimensional Langevin systems. The work is base on the theory presented in [VandenEijndenTransition2006, VandenEijndenTransition2010, VandenEijndenTowards2006](@cite) and [this python github repo](https://github.com/mar1akc/transition_path_theory_FEM_distmesh).

The functionality provided here is considered experimental and is currently restricted to two-dimensional Langevin dynamics of a particle exposed to a potential ``U``. Hence, the current API is not exported and should be loaded explicitly. We recommend checking out the [Transition Path Theory example](@ref TPT_example) for a demonstration of the available functionality.

```math
\begin{aligned}
\dot x &= p \,, \\
\dot p &= - \gamma p - \nabla U(x) + \sqrt{2\gamma\beta^{-1}} \dot W_t \,,
\end{aligned}
```
where ``x`` and ``p`` are the position and momentum of the particle, ``\gamma`` is the damping constant, ``\beta = 1/k_BT`` is the inverse temperature and ``W_t`` denotes a standard Wiener process.

To make the code more general, we need to implement the necessary methods for higher-dimensional systems and adapt the API accordingly. 

## Theory overview

TPT characterizes transition pathways between two metastable states A and B by computing:

- Forward and backward committor functions
- Reactive probability density
- Reactive current
- Transition rates

The calculations are performed on a triangulated mesh representing the system's state space.

```@docs
CriticalTransitions.LangevinSystem
```

### Committor functions

```@docs
CriticalTransitions.committor
```

### Invariant probability density

```@docs
CriticalTransitions.invariant_pdf
```

### Reactive current

```docs
CriticalTransitions.reactive_current
```

### Probability calculations

```@docs
CriticalTransitions.probability_reactive
```

```@docs
CriticalTransitions.probability_last_A
```

## References

```@bibliography
Pages = ["transition_path_theory.md"]
```