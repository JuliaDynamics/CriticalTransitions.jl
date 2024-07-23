# Quickstart

## Installation
As this module is not published yet, there are two ways to access it:

* **Option 1 (recommended): Install from GitHub**
    1. Enter the Julia package manager by typing `]` in the REPL: `julia> ]`
    2. type `add https://github.com/juliadynamics/CriticalTransitions.jl.git`
* **Option 2: Load module locally**
    1. Clone the repo: `git clone https://github.com/juliadynamics/CriticalTransitions.jl.git`
    2. In Julia, include the module file: `include("PATH/src/CriticalTransitions.jl")`, where `PATH` is the relative path to the repo you just cloned
    3. Load the module: `using .CriticalTransitions`

## Basic usage
The general workflow of `CriticalTransitions` essentially follows two steps:

1. Define your system (see [Define a CoupledSDEs system](@ref))
2. Investigate the system by calling methods (see [Methods](@ref))

!!! info "New system type: RateSystem"
    We are planning to introduce the the struct `RateSystem` along `CoupledSDEs`. In a `RateSystem`, the time dependence of parameters can conveniently be specified, laying the foundation for a toolbox to study rate-induced tipping, or R-tipping.

## Methods

Currently the following functions are implemented to analyze a [`CoupledSDEs`](@ref) and 
corresponding sample transition paths.

```@index
Pages = ["man/systemanalysis.md", "man/simulation.md", "man/sampling.md", "man/largedeviations.md"]
```