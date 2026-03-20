abstract type GMAMOptimizer end

"""
    GeometricGradient(; stepsize=0.01)

Projected-gradient descent optimizer for the geometric minimum action method
[heymann_pathways_2008](@citet).
"""
struct GeometricGradient{T} <: GMAMOptimizer
    stepsize::T
end
GeometricGradient(; stepsize=0.01) = GeometricGradient(stepsize)
