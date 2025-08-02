# Transition Path Theory

Transition Path Theory (TPT) provides a framework for analyzing rare transition events between metastable states in complex systems. This module implements TPT calculations for two-dimensional Langevin systems. The work is base on the theory presented in [vanden-eijnden_transition_2006, vanden-eijnden_towards_2006, vanden-eijnden_transition-path_2010](@cite).
The code is considered experimental as it only supports two-dimensional Langevin systems. To make the code more general, we need to implement the necessary methods for higher-dimensional systems and adapt the API accordingly.

## Theory Overview

TPT characterizes transition pathways between two metastable states A and B by computing:

- Forward and backward committor functions
- Reactive probability density
- Reactive current
- Transition rates

The calculations are performed on a triangulated mesh representing the system's state space.

```@docs
Langevin
```

### Committor Functions

```@docs
committor
```

### Invariant Probability Density

```@docs
invariant_pdf
```

### Reactive Current

```docs
reactive_current
```

### Probability Calculations

```@docs
probability_reactive
```

```@docs
probability_last_A
```

## References

```@bibliography
Pages = ["transition_path_theory.md"]
Canonical = false
```