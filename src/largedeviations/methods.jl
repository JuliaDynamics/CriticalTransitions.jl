abstract type GMAMOptimizer end

"""
$(TYPEDEF)

Projected-gradient descent optimizer for the geometric minimum action method
[heymann_pathways_2008](@citet).

# Fields
$(TYPEDFIELDS)

# Keyword constructors
$(METHODLIST)
"""
struct GeometricGradient{T <: Real} <: GMAMOptimizer
    """Step size for the projected gradient update."""
    stepsize::T
end

function GeometricGradient(; stepsize::Real = 0.01)
    T = typeof(float(stepsize))
    return GeometricGradient{T}(T(stepsize))
end
