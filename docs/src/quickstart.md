# Quickstart

## Installation
Although the package is not yet registered, you can install it from Github via the Julia package manager:
```julia
using Pkg; Pkg.add("https://github.com/juliadynamics/CriticalTransitions.jl.git")
```

The package is currently tested to be compatible with Julia versions `1.10` and `1.11`.

## Basic usage
The general workflow of `CriticalTransitions` essentially follows two steps:

1. Define your system as a `CoupledSDEs` (see [Define a CoupledSDEs system](@ref))
2. Investigate the system by calling methods (see [Methods](@ref))

!!! info "New system type: RateSystem"
    We are planning to introduce the the struct `RateSystem` along `CoupledSDEs`. In a `RateSystem`, the time dependence of parameters can conveniently be specified, laying the foundation for a toolbox to study rate-induced tipping, or R-tipping.

## Methods

Currently the following functions are implemented to analyze a [`CoupledSDEs`](@ref) and 
corresponding sample transition paths.

```@index
Pages = ["man/systemanalysis.md", "man/simulation.md", "man/sampling.md", "man/largedeviations.md", "man/utils.md"]
```