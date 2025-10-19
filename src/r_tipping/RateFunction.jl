"""
    RateFunction(;pfunc::Function, ptspan)

Time-dependent forcing protocol ``p(t)`` describing the evolution of a parameter over a
time interval `ptspan = (start, end)`. Used to construct a non-autonomous `RateSystem`.

## Keyword Arguments
- `pfunc`: function ``p(t)`` from ``ℝ → ℝ`` describing the parameter time dependence
- `ptspan`: domain interval over which `pfunc` is considered

## Description
The `RateFunction` type allows to specify the functional form of a parametric
forcing over a time interval. This forcing protocol can then be applied to the parameter
of a dynamical system using the `RateSystem` constructor.
The rate and magnitude of the forcing can be adjusted when constructing the `RateSystem`.

## Convenience functions
- `RateFunction(:linear)`: Creates a linear ramp from 0 to 1.
- `RateFunction(:tanh)`: Creates a hyperbolic tangent ramping from 0 to 1.

"""
struct RateFunction{F,T<:Real}
    pfunc::F
    ptspan::Tuple{T,T}
end

# Convenience constructors for RateFunction
function RateFunction(sym::Symbol)
    if sym == :linear
        return RateFunction(x -> x, (0.0, 1.0))
    elseif sym == :tanh
        return RateFunction(x -> tanh(x)/2 + 0.5, (-3.0, 3.0))
    else
        error("Only :linear or :tanh are supported input arguments.")
    end
end

# """
#     show(rc::RateFunction, n=50)

# Returns a `RateFunction` forcing protocol in the form (domain_values, pfunc_values)`,
# where `domain_values` are time values and `pfunc_values` the corresponding values of the
# forcing function `rc.pfunc`.

# ## Keyword arguments
# `n=50`: Number of data points to be returned.
# """
# function Base.show(rc::RateFunction, n::Integer=50)
#     domain_values = range(rc.interval[1], rc.interval[2]; length=n)
#     pfunc_values = rc.pfunc.(domain_values)  # broadcast
#     return domain_values, pfunc_values
# end
# TODO: use unicode plot to show how the rate function looks like
# https://docs.juliadsp.org/stable/windows/
