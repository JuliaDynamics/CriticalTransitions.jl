"""
    ForcingProfile(profile::Function, interval)

Time-dependent forcing profile ``p(t)`` describing the evolution of a parameter over a
domain interval `interval = (start, end)`. Used to define a parametric forcing when 
constructing a non-autonomous [`RateSystem`](@ref).

## Keyword arguments
- `profile`: function ``p(t)`` from ``ℝ → ℝ`` describing the parameter change
- `interval`: domain interval over which the `profile` is considered

Note: The units of `t` are arbitrary; the forcing profile is rescaled in system time units
and magnitude when using it to construct a `RateSystem`.

## Convenience functions
- `ForcingProfile(:linear)`: Creates a linear ramp from 0 to 1.
- `ForcingProfile(:tanh)`: Creates a hyperbolic tangent ramping from 0 to 1.
"""
struct ForcingProfile{F,T<:Real}
    profile::F
    interval::Tuple{T,T}
end

# Convenience constructors for ForcingProfile
function ForcingProfile(sym::Symbol)
    if sym == :linear
        return ForcingProfile(x -> x, (0.0, 1.0))
    elseif sym == :tanh
        return ForcingProfile(x -> tanh(x)/2 + 0.5, (-3.0, 3.0))
    else
        error("Only :linear or :tanh are supported input arguments.")
    end
end

# """
#     show(rc::ForcingProfile, n=50)

# Returns a `ForcingProfile` forcing protocol in the form (domain_values, pfunc_values)`,
# where `domain_values` are time values and `pfunc_values` the corresponding values of the
# forcing function `rc.profile`.

# ## Keyword arguments
# `n=50`: Number of data points to be returned.
# """
# function Base.show(rc::ForcingProfile, n::Integer=50)
#     domain_values = range(rc.interval[1], rc.interval[2]; length=n)
#     pfunc_values = rc.profile.(domain_values)  # broadcast
#     return domain_values, pfunc_values
# end
# TODO: use unicode plot to show how the rate function looks like
# https://docs.juliadsp.org/stable/windows/
