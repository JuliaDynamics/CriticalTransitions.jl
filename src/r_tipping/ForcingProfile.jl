"""
    ForcingProfile(profile::Function, interval)

Time-dependent forcing profile ``p(t)`` describing the evolution of a parameter over a
domain `interval = (start, end)`. Used to define a parametric forcing when
constructing a non-autonomous [`RateSystem`](@ref).

## Arguments

- `profile`: function ``p(t)`` from ``ℝ → ℝ`` describing the parameter change
- `interval`: domain interval `(start, end)` over which the `profile` is considered

Note: The units of ``t`` are arbitrary; the forcing profile is rescaled to system time units.

## Convenience functions

- `ForcingProfile(:linear)`: Create a linear ramp from 0 to 1.
- `ForcingProfile(:tanh)`: Create a hyperbolic tangent ramping
  from 0 to 1 with interval (-3, 3).
- [`data`](@ref): Get data points describing a given `ForcingProfile`.
"""
struct ForcingProfile{F, T <: Real}
    profile::F
    interval::Tuple{T, T}
end
ForcingProfile(x, y) = ForcingProfile(x, promote(y...))


# Convenience constructors for ForcingProfile
function ForcingProfile(sym::Symbol)
    if sym == :linear
        return ForcingProfile(x -> x, (0.0, 1.0))
    elseif sym == :tanh
        return ForcingProfile(x -> tanh(x) / 2 + 0.5, (-3.0, 3.0))
    else
        error("Only :linear or :tanh are supported input arguments.")
    end
end
