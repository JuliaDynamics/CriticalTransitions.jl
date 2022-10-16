# Getting started

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

## Methods

Currently the following functions are implemented to analyze a StochSystem.

```@index
Pages = ["man/simulation.md", "man/systemanalysis.md"]
```