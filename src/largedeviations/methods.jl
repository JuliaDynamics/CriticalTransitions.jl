abstract type GMAMOptimizer end

"""
$(TYPEDEF)

Optimizer configuration for the (s)gMAM projected-gradient update with built-in
backtracking line search (step size is provided per-call via keyword arguments).

# Fields
$(TYPEDFIELDS)

# Keyword constructors
$(METHODLIST)
"""
struct GeometricGradient{T<:Real} <: GMAMOptimizer
    """Step-size shrink factor on rejected steps (backtracking)."""
    shrink::T
    """Step-size growth factor after accepted steps (backtracking)."""
    grow::T
    """Maximum number of backtracking attempts per iteration."""
    max_backtracks::Int
    """Lower clamp for step size (backtracking)."""
    stepsize_min::T
    """Upper clamp for step size (backtracking)."""
    stepsize_max::T
end

function GeometricGradient(;
    shrink::Real=0.5,
    grow::Real=1.1,
    max_backtracks::Int=10,
    stepsize_min::Real=1e-12,
    stepsize_max::Real=1e3,
)
    T = promote_type(
        typeof(shrink),
        typeof(grow),
        typeof(stepsize_min),
        typeof(stepsize_max),
    )
    return GeometricGradient{T}(
        T(shrink),
        T(grow),
        max_backtracks,
        T(stepsize_min),
        T(stepsize_max),
    )
end
