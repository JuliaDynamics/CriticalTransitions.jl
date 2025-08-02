# Getting started

## Installation
To install the Julia language, we recommend [juliaup](https://github.com/JuliaLang/juliaup).

The `CriticalTransitions.jl` package can be installed from Github via the Julia package manager:
```julia
using Pkg; Pkg.add(url="https://github.com/juliadynamics/CriticalTransitions.jl.git")
```
You can then load the package with `using CriticalTransitions`.

The package is currently tested to be compatible with Julia versions `1.10` and `1.11`.

## Basic usage
The general workflow of `CriticalTransitions.jl` essentially follows two steps, similar to `DynamicalSystems.jl`:

1. Define your dynamical system (see [Defining a forced dynamical system](@ref))
2. Investigate the system by calling functions (see [Index](@ref))

Some functions are only loaded as extensions when loading other dependency packages (see [Extensions](@ref)).

## Documentation
The [Tutorial](@ref) and code examples in the *Examples* section illustrate some use cases of the package. All available functions and types are documented in the *Manual* section (see [Index](@ref) for an overview).

## Index
**Defining a system and its forcing**
- [`DynamicalSystemsBase.CoupledODEs`](@ref)
- [`DynamicalSystemsBase.CoupledSDEs`](@ref) (alias for `StochSystem`)

**System analysis and simulation**
```@index
Pages = ["man/systemanalysis.md", "man/simulation.md"]
```

**Stochastic dynamics: Transition path sampling, large deviation theory**
```@index
Pages = ["man/sampling.md", "man/largedeviations.md"]
```

**Nonautonomous dynamics: R-tipping** 
```@index
Pages = ["man/r-tipping.md"]
```

**More**
```@index
Pages = ["man/utils.md"]
```

## Extensions
- Loading the `ChaosTools` dependency makes the functions [`fixedpoints`](@ref) and [`intervals_to_box`](@ref) available.

```julia
using CriticalTransitions
using ChaosTools

# Now the following functions are available
fixedpoints()
intervals_to_box()
```

- Loading the `Attractors` dependency makes the functions [`edgetracking`](@ref) and [`bisect_to_edge`](@ref) available.

```julia
using CriticalTransitions
using Attractors

# Now the following functions are available
edgetracking()
bisect_to_edge()
```