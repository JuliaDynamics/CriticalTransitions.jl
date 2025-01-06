# Transition Path Theory

Transition Path Theory (TPT) provides a framework for analyzing rare transition events between metastable states in complex systems. This module implements TPT calculations for Langevin dynamics systems.

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

1. W. E and E. Vanden-Eijnden, "Towards a Theory of Transition Paths", Journal of Statistical Physics, 2006
