"""
    RateConfig(pfunc::Function, interval)

Time-dependent forcing protocol ``p(t)`` describing the evolution of a parameter over a
time interval `interval = (start, end)`.
Used to construct a non-autonomous `RateSystem`.

## Arguments
- `pfunc`: function ``p(t)`` from ``ℝ → ℝ`` describing the parameter time dependence
- `interval`: domain interval over which `pfunc` is considered

## Description
The `RateConfig` type allows to specify the functional form of a parametric
forcing over a time interval. This forcing protocol can then be applied to the parameter
of a dynamical system using the `RateSystem` constructor, which also allows to modify
the rate and magnitude of the forcing.

## Convenience functions
`RateConfig(:linear)`: Creates a linear ramp from 0 to 1.
`RateConfig(:tanh)`: Creates a hyperbolic tangent ramping from 0 to 1.

Use `show(::RateConfig)` to display ``t`` and ``p(t)``.
"""
struct RateConfig{F,T<:Real}
    pfunc::F
    interval::Tuple{T,T}
end

# Convenience constructors for RateConfig
function RateConfig(sym::Symbol)
    if sym == :linear
        return RateConfig(x -> x, (0.0, 1.0))
    elseif sym == :tanh
        return RateConfig(x -> tanh(x)/2 + 0.5, (-3.0, 3.0))
    else
        error("Only :linear or :tanh are supported input arguments.")
    end
end

"""
    show(rc::RateConfig, n=50)

Returns a `RateConfig` forcing protocol in the form (domain_values, pfunc_values)`,
where `domain_values` are time values and `pfunc_values` the corresponding values of the
forcing function `rc.pfunc`.

## Keyword arguments
`n=50`: Number of data points to be returned.
"""
function show(rc::RateConfig, n::Integer=50)
    domain_values = range(rc.interval[1], rc.interval[2]; length=n)
    pfunc_values = rc.pfunc.(domain_values)  # broadcast
    return domain_values, pfunc_values
end
