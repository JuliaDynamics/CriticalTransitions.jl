# TODO: Generalize to multi-parameter
# In general, all instances of parameters that are just `Real`
# must be made into `Vector{Pair{IDX, Real}}` with `IDX` the index type
# This also decreases the number of fields as we dont need different fields
# for different parameters. we can have a single `pspec` for pawrameter specification
# that is `Vector{Pair{IDX, Real}}`. Or it can be a `Vector{Dict}`. In
# general anything that can be given to `set_parameters!`. So we cross ref that
# docstring and everything becomes simple and unified and harmonious...

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
    #"Parameter index of the time-dependent parameter"
    #pidx::Vector{Pair{P, T}}
    "Forcing profile"
    fp::Vector{ForcingProfile{F,T}}
    "Start time of parametric forcing"
    forcing_start_time::T
    "Duration of parametric forcing"
    forcing_duration::T
    "Magnitude of parametric forcing"
    forcing_scale::T
    "Parameter indices and their initial (autonomous) value(s)"
    p0::Vector{Pair{P, T}} # make this a vector
    "Initial time (of system initiation)"
    t0::T
end

"""
    RateSystem <: ContinuousTimeDynamicalSystem
    RateSystem(ds::ContinuousTimeDynamicalSystem, fp::ForcingProfile, pidx; kwargs...)

A non-autonomous dynamical system constructed by applying a time-dependent parametric
forcing protocol ([`ForcingProfile`](@ref)) to an underlying autonomous continuous-time
dynamical system `ds`.

The `ForcingProfile`, applied to the parameter with index `pidx`, describes the functional
form of the parameter evolution over a given `interval`. This interval is rescaled
in system time units by defining its start time (`forcing_start_time`) and duration
(`forcing_duration`). The magnitude of the forcing can be adjusted by scaling the
forcing profile by a factor `forcing_scale`.

The integer `pidx` refers to the relevant index of the parameter container `p` of `ds`.

Before `t = forcing_start_time`, the system parameters correspond to that of the
underlying autonomous system `ds`. At the end of the forcing interval
(`t = forcing_start_time + forcing_duration`), the parameters are frozen at their current
value and thereafter the system is again autonomous.

## Keyword arguments
- `forcing_start_time = fp.interval[1]`: Time when parameter change starts
- `forcing_duration = fp.interval[2] - fp.interval[1]`: Duration of the parameter change (in system time units)
- `forcing_scale = 1.0`: Amplitude multiplication factor of the forcing protocol
- `t0 = 0.0`: Initial time of the system
"""
struct RateSystem{S,F,P,T} <: ContinuousTimeDynamicalSystem
    "Non-autonomous continuous-time dynamical system"
    system::S
    "Data structure representing the system and forcing specifications in system units"
    forcing::RateSystemSpecs{F,P,T}
end

function RateSystem(
    ds::ContinuousTimeDynamicalSystem,
    fp::Vector{ForcingProfile},
    pidx::Vector{P};
    forcing_start_time::T=initial_time(ds),
    forcing_duration::T=(fp.interval[2] - fp.interval[1]),
    forcing_scale::T=1.0,
    t0::T=initial_time(ds),
) where {P,T<:Real}
    (forcing_start_time >= t0) ||
        throw("The forcing cannot start before the system's initial time `t0`
              but your forcing_start_time ($(forcing_start_time)) < t0 ($(t0)).")

    # TODO: Make p0 a vector
    p0 = [k => (current_parameter(ds, k)) for k in 1:length(pidx)]
    rss = RateSystemSpecs(
        dynamic_rule(ds),
        fp,
        forcing_start_time,
        forcing_duration,
        forcing_scale,
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

function (rss::RateSystemSpecs)(u, p, t)
    # TODO: this should be rewritten with `current_parameter`
    p_at_t = p_modified(rss, t)
    for k in 1:length(rss.p0)
        pidx, pval = p_at_t[k]
        if p isa Union{AbstractArray,AbstractDict}
            setindex!(p, pval, pidx)
        else
            setproperty!(p, pidx, pval)
        end
    end
    # TODO: This doesn't work if the dynamical system is in-place...
    # We need to multi-dispatch one more version of `(rss::)` function with 4 arguments
    return rss.unforced_rule(u, p, rss.t0)
end

# Returns the non-autonomous system's parameter value at time t
function p_modified(rss::RateSystemSpecs, t::Real)
    # TODO: Make `p0` vector, and make the function re-write a dummy vector of time dependent parameters
    p_timedep = Vector{typeof(rss.p0)}(undef, length(rss.p0))

    for k in 1:length(rss.p0)
        pidx, p0 = rss.p0[k]
        f = rss.fp[k].profile
        fstart = rss.fp[k].interval[1]
        fend = rss.fp[k].interval[2]
        if t <= rss.forcing_start_time
            p_timedep[k] = pidx => p0
        else
            if t < rss.forcing_start_time + rss.forcing_duration
                time_shift = ((fend - fstart)/rss.forcing_duration)*(t -
                    rss.forcing_start_time) + fstart
                p_timedep[k] = pidx => p0 + rss.forcing_scale*(f(time_shift) - f(fstart))
            else
                p_timedep[k] = pidx => p0 + rss.forcing_scale*(f(fend) - f(fstart))
            end
        end
    end
    return p_timedep
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
$(TYPEDSIGNATURES)

Returns the `i`-th parameter of a [`RateSystem`](@ref) at time `t` (in system time units).
"""
function parameter(rs::RateSystem, i, t)
    p = deepcopy(current_parameter(rs, i))
    p_mod = p_modified(rs.forcing, t)
    timedependent_indices = [pmod[k][1] for k in 1:length(pmod)]
    (i in timedependent_indices) ? (return p_mod[i][2]) : (return p[i])
end

"""
$(TYPEDSIGNATURES)

Returns the parameter vector of a [`RateSystem`](@ref) at time `t` (in system time units).
"""
function parameters(rs::RateSystem, t)
    p = deepcopy(current_parameters(rs))
    for (pidx, pval) in p_modified(rs.forcing, t)
        p[pidx] = pval
    end
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
