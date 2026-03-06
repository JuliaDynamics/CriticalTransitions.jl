"""
    ForcingProfile(profile::Function, interval)

Time-dependent forcing profile ``p(t)`` describing the evolution of a parameter over a
domain interval `interval = (start, end)`. Used to define a parametric forcing when
constructing a non-autonomous [`RateSystem`](@ref).

## Keyword arguments

- `profile`: function ``p(t)`` from ``ℝ → ℝ`` describing the parameter change
- `interval`: domain interval over which the `profile` is considered

Note: The units of ``t`` are arbitrary; the forcing profile is rescaled to system time units.

## Convenience functions

- `ForcingProfile(:linear)`: Create a linear ramp from 0 to 1.
- `ForcingProfile(:tanh)`: Create a hyperbolic tangent ramping
  from 0 to 1 with interval (-3, 3).
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

"""
    RateSystemSpecs <: Function

A mutable data structure storing information needed to construct and modify a
[`RateSystem`](@ref).

Calling `(::RateSystemSpecs)(u,p,t)` returns the nonautonomous drift of the `RateSystem`
at time `t`, where `p` is the parameter container of the underlying autonomous system.
"""
mutable struct RateSystemSpecs{R,F,P,T} <: Function
    "Dynamic rule of the underlying autonomous system"
    unforced_rule::R
    "Parameter index of the time-dependent parameter"
    pidx::P # make this a vector
    "Forcing profile"
    fp::ForcingProfile{F,T}
    "Start time of parametric forcing"
    forcing_start_time::T
    "Duration of parametric forcing"
    forcing_duration::T
    "Magnitude of parametric forcing"
    forcing_scale::T
    "Initial (autonomous) parameter value"
    p0::T # make this a vector
    "Initial time (of system initiation)"
    t0::T
end

"""
    RateSystem <: ContinuousTimeDynamicalSystem
    RateSystem(ds::ContinuousTimeDynamicalSystem, fp::ForcingProfile, pidx; kw...)

A non-autonomous dynamical system constructed by applying a time-dependent parametric
forcing protocol(s) to an underlying autonomous continuous-time dynamical system `ds`.

The parameter of the system with index `pidx` will be forced according to the given
[`ForcingProfile`](@ref). The `ForcingProfile` describes the functional form of the parameter
evolution over a given `interval`.

## Keyword arguments

- `forcing_start_time = fp.interval[1]`: Time when parameter change starts.
- `forcing_duration = fp.interval[2] - fp.interval[1]`:
   Duration of the parameter change (in system time units).
- `forcing_scale = 1.0`: Amplitude multiplication factor of the forcing protocol.
- `t0 = 0.0`: Initial time of the system.

## Description

The time interval (in units of `ds`) where forcing applied is from
`forcing_start_time` to `forcing_duration`. As such, the `interval` given
of the [`ForcingProfile`](@ref) is rescaled in system time units; its purpose is only
to establish the functional form of the forcing.
The parameter value starts at `p0 = current_parameter(ds, pidx)`.
Its maximum value is given by the maximum of the forcing function given to the
[`ForcingProfile`](@ref), multiplied by `forcing_scale`.

Before `t = forcing_start_time`, the system parameters correspond to that of the
underlying autonomous system `ds`. At the end of the forcing interval
(`t = forcing_start_time + forcing_duration`), the parameters are frozen at their final
value and thereafter the system is again autonomous.

## Multiple parameters

    RateSystem(ds::ContinuousTimeDynamicalSystem, forcing_profiles::Dict; kw...)

Use the above signature with `forcing_profiles` a dictionary mapping parameter indices
(anything valid for `set_parameter`) to [`ForcingProfile`](@ref) instances.

    RateSystem(ds::ContinuousTimeDynamicalSystem, fp::ForcingProfile, pidxs::Vector; kw...)

Use this when all parameters share the same forcing profile and pass all parameter
indices to `pidxs`.
"""
struct RateSystem{S,F,P,T} <: ContinuousTimeDynamicalSystem
    "Non-autonomous continuous-time dynamical system"
    system::S
    "Data structure representing the system and forcing specifications in system units"
    forcing::RateSystemSpecs{F,P,T}
end

function RateSystem(
        ds::ContinuousTimeDynamicalSystem,
        fp::ForcingProfile,
        pidx;
        forcing_start_time=initial_time(ds),
        forcing_duration=(fp.interval[2] - fp.interval[1]),
        forcing_scale=1.0,
        t0=initial_time(ds),
    )
    (forcing_start_time >= t0) ||
        throw("The forcing cannot start before the system's initial time `t0`
              but your forcing_start_time ($(forcing_start_time)) < t0 ($(t0)).")

    # TODO: allow multi parameter
    p0 = current_parameter(ds, pidx)
    # promote
    a, b, c = float.((forcing_start_time, forcing_duration, forcing_scale))
    rss = RateSystemSpecs(
        dynamic_rule(ds),
        pidx,
        fp,
        a,
        b,
        c,
        p0,
        t0,
    )

    # TODO: Do we want the rate system to have the same parameter container,
    # or a deep copy of it...? Derivative systems in DynamicalSystems.jl
    # typically have the same parameter container, however this one is special,
    # as it constantly alters the parameter container...
    newds = CoupledODEs(rss, current_state(ds), deepcopy(current_parameters(ds)); t0)
    return RateSystem(newds, rss)
end

# TODO: this must be rewritten using `set_parameter!` or its source code.
# otherwise it doens't work with ModelingToolkit.jl
function (rss::RateSystemSpecs)(u, p, t)
    p_at_t = p_modified(rss, t)
    if p isa Union{AbstractArray,AbstractDict}
        setindex!(p, p_at_t, rss.pidx)
    else
        setproperty!(p, rss.pidx, p_at_t)
    end
    # TODO: This doesn't work if the dynamical system is in-place...
    # We need to multi-dispatch one more version of `(rss::)` function with 4 arguments
    return rss.unforced_rule(u, p, rss.t0)
end

# Returns the non-autonomous system's parameter value at time t
function p_modified(rss::RateSystemSpecs, t::Real)
    # TODO: Make `p0` vector, and make the function re-write a dummy vector of time dependent parameters
    p0 = rss.p0
    f = rss.fp.profile
    section_start = rss.fp.interval[1]
    section_end = rss.fp.interval[2]
    # making the function piecewise constant with range [p0, p0+forcing_scale]
    if t ≤ rss.forcing_start_time
        return p0
    else
        if t < rss.forcing_start_time + rss.forcing_duration
            # performing the time shift corresponding to stretching/squeezing
            time_shift =
                ((section_end - section_start) / rss.forcing_duration) *
                (t - rss.forcing_start_time) + section_start
            return p0 + rss.forcing_scale*(f(time_shift) - f(section_start))
        else
            return p0 + rss.forcing_scale*(f(section_end) - f(section_start))
        end
    end
end

"""
$(TYPEDSIGNATURES)

Sets the amplitude (`forcing_scale`) of the forcing protocol applied to
the [`RateSystem`](@ref) `rs`.
"""
function set_forcing_scale!(rs::RateSystem, scale)
    rs.forcing.forcing_scale = scale
    return rs
end

"""
$(TYPEDSIGNATURES)

Sets the duration (`forcing_duration`) of the forcing protocol applied to
the [`RateSystem`](@ref) `rs`.
"""
function set_forcing_duration!(rs::RateSystem, duration)
    rs.forcing.forcing_duration = duration
    return rs
end

"""
$(TYPEDSIGNATURES)

Sets the start time (`forcing_start_time`) of the forcing protocol applied to
the [`RateSystem`](@ref) `rs`.
"""
function set_forcing_start!(rs::RateSystem, start_time)
    rs.forcing.forcing_start = start_time
    return rs
end

"""
    current_parameters(rs::RateSystem, t)

Returns the parameter vector of a [`RateSystem`](@ref) at time `t` (in system time units).
Note: this function returns a copy of the original parameter container.
"""
function DynamicalSystemsBase.current_parameters(rs::RateSystem, t)
    p = deepcopy(current_parameters(rs))
    # TODO: Doesn't work for struct parameters, use generic function from DynamicalSystemsBase
    p[rs.forcing.pidx] = p_modified(rs.forcing, t)
    return p
end

"""
$(TYPEDSIGNATURES)

Returns an autonomous dynamical system of type [`DynamicalSystemsBase.CoupledODEs`](@ref)
corresponding to the frozen system of the non-autonomous [`RateSystem`](@ref) `rs` at
time `t`.
"""
function frozen_system(rs::RateSystem, t)
    ds = CoupledODEs(rs.forcing.unforced_rule, current_state(rs), parameters(rs, t))
    return ds
end

# Extensions
for f in (
    :initial_state,
    :initial_parameters,
    :current_parameter,
    :current_parameters,
    :dynamic_rule,
    :current_time,
    :current_state,
    :initial_time,
    :successful_step,
    :set_parameter!,
    :set_parameters!,
    :trajectory,
) # all api functions here
    @eval DynamicalSystemsBase.$(f)(rs::RateSystem, args...; kw...) = $(f)(
        rs.system, args...; kw...
    )
end

reinit!(rs::RateSystem, u=initial_state(rs); kw...) = reinit!(rs.system, u; kw...)

function DynamicalSystemsBase.set_state!(rs::RateSystem, u::AbstractArray)
    return (DynamicalSystemsBase.set_state!(rs.system, u); rs)
end

SciMLBase.step!(rs::RateSystem, args...) = (SciMLBase.step!(rs.system, args...); rs)
SciMLBase.isinplace(rs::RateSystem) = SciMLBase.isinplace(rs.system)
StateSpaceSets.dimension(rs::RateSystem) = StateSpaceSets.dimension(rs.system)
DynamicalSystemsBase.isdeterministic(rs::RateSystem) = true # TODO: only true if CoupledODEs

(rs::RateSystem)(t::Real) = rs.system(t)
