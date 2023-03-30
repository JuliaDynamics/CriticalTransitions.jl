# CriticalTransitions.jl

[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://reykboerner.github.io/CriticalTransitions.jl/dev/) [![](https://img.shields.io/badge/docs-stable-green.svg)](https://reykboerner.github.io/CriticalTransitions.jl/stable/)

A Julia package for the numerical investigation of **noise- and rate-induced transitions in dynamical systems**.

Building on [DynamicalSystems.jl](https://juliadynamics.github.io/DynamicalSystems.jl/stable/) and [DifferentialEquations.jl](https://diffeq.sciml.ai/stable/), this package aims to provide a toolbox for dynamical systems under time-dependent forcing, with a focus on tipping phenomena and metastability.
## Usage
See [package documentation](https://reykboerner.github.io/CriticalTransitions.jl/stable/).

## Example: Bistable FitzHugh-Nagumo model
```julia
using CriticalTransitions

# System parameters
p = [1., 3., 1., 1., 1., 0.]
noise_strength = 0.02

# Define stochastic system
sys = StochSystem(fitzhugh_nagumo, p, zeros(2), noise_strength)

# Get stable fixed points
fps, eigs, stab = fixedpoints(sys, [-2,-2], [2,2])
fp1, fp2 = fps[stab]

# Generate noise-induced transition from one fixed point to the other
path, times, success = transition(sys, fp1, fp2)

# ... and more, check out the documentation!
```

---

Developers: Reyk Börner, Ryan Deeley and Raphael Römer

Thanks to Jeroen Wouters, Calvin Nesbitt, Tobias Grafke, George Datseris and Oliver Mehling

This work is part of the [CriticalEarth](https://criticalearth.eu) project.