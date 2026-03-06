# # CriticalTransitions.jl Tutorial

# !!! note "DynamicalSystems.jl and Attractors.jl background recommended"
#     CriticalTransitions.jl is an advanced software for the analysis of critical trasitions
#     in dynamical systems. Due to its advanced nature it is recommended that you have basic
#     familiarity with the DynamicalSystems.jl and Attractors.jl packages, by going through
#     their main tutorials.

# The general workflow of CriticalTransitions.jl consists of two steps, similar to DynamicalSystems.jl:

# 1. Define your special dynamical system type
#    (either [`RateSystem`](@ref) or [`StochasticSystem`](@ref), see below in this tutorial).
# 2. Investigate the system by calling existing functions on it (see [API](@ref) and examples in this tutorial).

# The picture below showcases the two main routes one can go: rate or stochastic



# ## Creating dynamical systems for critical transitions



