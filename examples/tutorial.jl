# # CriticalTransitions.jl Tutorial

# !!! note "DynamicalSystems.jl and Attractors.jl background recommended"
#     CriticalTransitions.jl is an advanced software for the analysis of critical trasitions
#     in dynamical systems. Due to its advanced nature it is recommended that you have basic
#     familiarity with the DynamicalSystems.jl and Attractors.jl packages, by going through
#     their main tutorials.

# The general workflow of CriticalTransitions.jl consists of two steps, similar to DynamicalSystems.jl:

# 1. Define your special dynamical system type
#    (either [`RateSystem`](@ref) or [`RandomSystem`](@ref), see below in this tutorial).
# 2. Investigate the system by calling existing functions on it (see [API](@ref) and examples in this tutorial).

# The picture below showcases the two main routes one can go: rate or stochastic

# ## Types of systems

# There are two main system types supported by this package, both being non-autonomous.
# The first type is systems that are driven by noise, primarily stochastic ordinary differential equations.
# The second type is systems whose parameters change with time deterministically according.

# In both cases one often starts with an autonomous deterministic system.
# For the scope of this tutorial, this will be the FitzHugh-Nagumo model:

# ```math
# \begin{aligned}
# \frac{du}{dt} &= \frac{1}{\epsilon} \left( -\alpha u^3 + \gamma u - \kappa v + I \right) \\
# \frac{dv}{dt} &= -\beta v + u \ ,
# \end{aligned}
# ```


# ## RateSystem: creation

# Transforming a deterministic `DynamicalSystem` to a `RateSystem` is straightforward.
# All we have to do is define a number of profiles

...

# ## RateSystem: example applications

# ## RandomSystem: creation

# ## RandomSystem: example applications
