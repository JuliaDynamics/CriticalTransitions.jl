# Quickstart

## Installation
As this module is not published yet, there are two ways to access it:

* **Install from GitHub**
    1. Enter the Julia package manager `julia> ]`
    2. type `add https://github.com/reykboerner/CriticalTransitions.jl.git`
* **Load module locally**
    1. Clone the repo: `git clone https://github.com/reykboerner/CriticalTransitions.jl.git`
    2. In Julia, include the module file: `include("PATH/src/CriticalTransitions.jl")`, where `PATH` is the relative path to the repo you just cloned
    3. Load the module: `using .CriticalTransitions`

## Basic usage
The general workflow of `CriticalTransitions` essentially follows two steps:

1. Define your system (see [Defining a StochSystem](@ref))
2. Investigate the system by calling methods (see [Methods](@ref))

!!! info "Extension to RateSystem and TippingSystem"
    We are currently working on extending the types of dynamical systems that can be studied with CriticalTransitions.jl. Particularly, we are planning to introduce the overarching structure `TippingSystem`, which has two subtypes: `StochSystem` (as it already exists) and `RateSystem`, a new dynamical system type in which the system parameters may evolve in time.

## Methods

Currently the following functions are implemented to analyze a [`StochSystem`](@ref) and 
corresponding sample transition paths.

```@index
Pages = ["man/systemanalysis.md", "man/simulation.md", "man/sampling.md", "man/largedeviations.md"]
```

## Systems

The following deterministic ODE systems have been implemented so far:

```@index
Pages = ["man/systems.md"]
```