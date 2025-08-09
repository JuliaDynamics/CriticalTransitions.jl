# Sampling transitions

## Direct simulation
These functions generate noise-induced transitions between a given initial and final state
via Monte Carlo rejection sampling. While [`transition`](@ref) samples one transition,
the [`transitions`](@ref) function generates an ensemble of transition paths. It supports
parallelization on multiple threads to speed up the run.

```@docs
transition
transitions
```

Results from the [`transitions`](@ref) function are organized in the following types:

```@docs
CriticalTransitions.TransitionEnsemble
CriticalTransitions.TransitionStatistics
```

## Pathspace sampling

!!! todo "System-independent implementation of pathspace samping"
    We are planning to add a general implementation of pathspace sampling as described
    in [borner2024saddle](@citet) (see article's Supplemental Material), where transition
    paths are sampled as stochastic bridges between a fixed start and end point in state
    space under stochastic gradient descent in pathspace.