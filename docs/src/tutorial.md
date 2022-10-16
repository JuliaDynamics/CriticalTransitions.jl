# Tutorial

## Example: Bistable FitzHugh-Nagumo model
```julia
include("src/CriticalTransitions.jl")
using .CriticalTransitions

# Parameters
p = [1., 3., 1., 1., 1., 0.]

# StochSystem
sys = StochSystem(FitzHughNagumo, p, 2)

# Get stable fixed points
eqs, eigs, stab = fixedpoints(sys, [-2,-2], [2,2])
fp1, fp2 = eqs[stab]

# Simulate noisy trajectory starting from fixed point 1
sim = simulate(sys, fp1)
```