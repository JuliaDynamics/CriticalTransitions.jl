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
mutable struct RateSystemSpecs{R,K,T,PC,E} <: Function
    "Dynamic rule of the underlying autonomous system"
    unforced_rule::R
    "Mapping parameter index => ForcingProfile"
    forcers::Dict{K,ForcingProfile}
    "Mapping parameter index => forcing start time"
    forcing_start_time::Dict{K,T}
    "Mapping parameter index => forcing duration"
    forcing_duration::Dict{K,T}
    "Mapping parameter index => forcing scale"
    forcing_scale::Dict{K,T}
    "Mapping parameter index => initial (autonomous) parameter value"
    p0::Dict{K,E}
    "Initial time (of system initiation)"
    t0::T
    "Dummy container (copy of parameters)"
    pdummy::PC
    "Reference to the owning (Coupled) system, set after construction"
    owner::Any
end

# Accessors for initial parameter(s). These are explicit and avoid overriding `getproperty`.
function initial_parameter_map(rss::RateSystemSpecs)
    return getfield(rss, :p0)
end

function initial_parameter(rss::RateSystemSpecs, pkey)
    return getfield(rss, :p0)[pkey]
end

function initial_parameter_value(rss::RateSystemSpecs)
    p0map = getfield(rss, :p0)
    if length(getfield(rss, :forcers)) == 1
        return first(values(p0map))
    else
        throw(ArgumentError("initial_parameter_value(rss) only valid when exactly one forcer is present"))
    end
end

# RateSystem-level wrappers


# `parameters` convenience wrapper (moved here from the package entry):
function parameters(sys, args...; kw...)
    return current_parameters(sys, args...; kw...)
end

"""
    RateSystem <: ContinuousTimeDynamicalSystem
    RateSystem(ds::ContinuousTimeDynamicalSystem, fp::ForcingProfile, pidx; kw...)

A non-autonomous dynamical system constructed by applying a time-dependent parametric
forcing protocol(s) to an underlying autonomous continuous-time dynamical system `ds`.

The parameter of the system with index `pidx` will be forced according to the given
[`ForcingProfile`](@ref). The `ForcingProfile` describes the functional form of the parameter
evolution over a given `interval`. See below for forcing multiple parameters.

## Keyword arguments

- `forcing_start_time = fp.interval[1]`: Time when parameter change starts.
- `forcing_duration = fp.interval[2] - fp.interval[1]`:
   Duration of the parameter change (in system time units).
- `forcing_scale = 1.0`: Amplitude multiplication factor of the forcing protocol.
- `t0 = 0.0`: Initial time of the system.

The functions [`set_forcing_start!`](@ref), [`set_forcing_duration`](@ref), and
[`set_forcing_scale`](@ref) can be used to update the forcing after system creation.

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
(anything valid for `set_parameter!`) to [`ForcingProfile`](@ref) instances, one
for each parameter. In this scenario the keywords `forcing_start_time, forcing_duration,
forcing_scale` also become dictionaries mapping parameter indices to the respective values.
"""
struct RateSystem{S,R} <: ContinuousTimeDynamicalSystem
    "Non-autonomous continuous-time dynamical system"
    system::S
    "Data structure representing the system and forcing specifications in system units"
    forcing::RateSystemSpecs{R}
end

# RateSystem-level wrappers (placed after `RateSystem` type definition)
function initial_parameter_map(rs::RateSystem)
    return initial_parameter_map(rs.forcing)
end

function initial_parameter(rs::RateSystem, pkey)
    return initial_parameter(rs.forcing, pkey)
end

function initial_parameter_value(rs::RateSystem)
    return initial_parameter_value(rs.forcing)
end

function RateSystem(
        ds::ContinuousTimeDynamicalSystem,
        forcer::Dict;
        forcing_start_time=nothing,
        forcing_duration=nothing,
        forcing_scale=nothing,
        t0=initial_time(ds),
    )

    # Normalize inputs into per-parameter dictionaries
    forcers = Dict(forcer) # shallow copy

    # start times
    start_map = Dict{Any,Real}()
    if forcing_start_time === nothing
        for (k, _) in forcers
            start_map[k] = t0
        end
    elseif isa(forcing_start_time, AbstractDict)
        start_map = Dict(forcing_start_time)
    else
        for (k, _) in forcers
            start_map[k] = float(forcing_start_time)
        end
    end

    # durations
    duration_map = Dict{Any,Real}()
    if forcing_duration === nothing
        for (k, fp) in forcers
            duration_map[k] = float(fp.interval[2] - fp.interval[1])
        end
    elseif isa(forcing_duration, AbstractDict)
        duration_map = Dict(forcing_duration)
    else
        for (k, _) in forcers
            duration_map[k] = float(forcing_duration)
        end
    end

    # scales
    scale_map = Dict{Any,Real}()
    if forcing_scale === nothing
        for (k, _) in forcers
            scale_map[k] = 1.0
        end
    elseif isa(forcing_scale, AbstractDict)
        scale_map = Dict(forcing_scale)
    else
        for (k, _) in forcers
            scale_map[k] = float(forcing_scale)
        end
    end

    # initial autonomous parameter values
    if isempty(forcers)
        throw(ArgumentError("`forcer` dictionary must contain at least one entry"))
    end

    p0_map = Dict{Any,Any}()
    for (k, _) in forcers
        p0_map[k] = current_parameter(ds, k)
    end

    # determine concrete parameter types for strong typing
    R = typeof(dynamic_rule(ds))
    first_key = first(keys(forcers))
    K = typeof(first_key)
    # infer time type from first ForcingProfile interval
    first_profile = first(values(forcers))
    Ttype = typeof(first_profile.interval[1])
    PC = typeof(current_parameters(ds))
    E = typeof(p0_map[first_key])

    pdummy = deepcopy(current_parameters(ds))

    rss = RateSystemSpecs{R,K,Ttype,PC,E}(
        dynamic_rule(ds),
        forcers,
        start_map,
        duration_map,
        scale_map,
        p0_map,
        t0,
        pdummy,
        nothing,
    )

    newds = CoupledODEs(rss, current_state(ds), deepcopy(current_parameters(ds)); t0)
    # set owner and ensure pdummy has same shape as system parameters
    rss.owner = newds
    rss.pdummy = deepcopy(current_parameters(newds))
    return RateSystem(newds, rss)
end

# TODO: this must be rewritten using `set_parameter!` or its source code.
# otherwise it doens't work with ModelingToolkit.jl;
# Or better yet, use `set_parameters!` and give a dict of parameters to set?
function (rss::RateSystemSpecs)(u, p, t)
    pmod = p_modified(rss, p, t)
    # Prefer calling out-of-place signature `f(u,p,t)` which returns du
    return rss.unforced_rule(u, pmod, t)
end

function (rss::RateSystemSpecs)(du, u, p, t)
    pmod = p_modified(rss, p, t)
    # Try calling in-place signature `f(du,u,p,t)` first; if not defined,
    # fall back to out-of-place `f(u,p,t)` and copy into `du`.
    try
        return rss.unforced_rule(du, u, pmod, t)
    catch err
        # Fall back to out-of-place: get result and copy
        out = rss.unforced_rule(u, pmod, t)
        try
            du .= out
            return nothing
        catch
            # If broadcasting assignment fails, rethrow original method error
            rethrow(err)
        end
    end
end

# Returns the non-autonomous system's parameter value at time t
function p_modified(rss::RateSystemSpecs, p, t::Real)
    # In-place update for Dict and Vector parameter containers (mutates `p`),
    # otherwise attempt to update the owning system via `set_parameter!`/`set_parameters!`.
    if isa(p, AbstractDict)
        for (pkey, profile) in rss.forcers
            p0 = getfield(rss, :p0)[pkey]
            f = profile.profile
            section_start = profile.interval[1]
            section_end = profile.interval[2]
            start_time = rss.forcing_start_time[pkey]
            duration = rss.forcing_duration[pkey]
            scale = rss.forcing_scale[pkey]

            if t <= start_time
                pt = p0
            elseif t < start_time + duration
                time_shift = ((section_end - section_start) / duration) * (t - start_time) + section_start
                pt = p0 + scale * (f(time_shift) - f(section_start))
            else
                pt = p0 + scale * (f(section_end) - f(section_start))
            end

            p[pkey] = pt
        end
        return p

    elseif isa(p, AbstractVector)
        for (pkey, profile) in rss.forcers
            idx = Int(pkey)
            p0 = getfield(rss, :p0)[pkey]
            f = profile.profile
            section_start = profile.interval[1]
            section_end = profile.interval[2]
            start_time = rss.forcing_start_time[pkey]
            duration = rss.forcing_duration[pkey]
            scale = rss.forcing_scale[pkey]

            if t <= start_time
                pt = p0
            elseif t < start_time + duration
                time_shift = ((section_end - section_start) / duration) * (t - start_time) + section_start
                pt = p0 + scale * (f(time_shift) - f(section_start))
            else
                pt = p0 + scale * (f(section_end) - f(section_start))
            end

            p[idx] = pt
        end
        return p
    else
        # Arbitrary parameter container: build a copy in pd and try to apply updates
        pd = deepcopy(rss.pdummy)
        for (pkey, profile) in rss.forcers
            p0 = getfield(rss, :p0)[pkey]
            f = profile.profile
            section_start = profile.interval[1]
            section_end = profile.interval[2]
            start_time = rss.forcing_start_time[pkey]
            duration = rss.forcing_duration[pkey]
            scale = rss.forcing_scale[pkey]

            if t <= start_time
                pt = p0
            elseif t < start_time + duration
                time_shift = ((section_end - section_start) / duration) * (t - start_time) + section_start
                pt = p0 + scale * (f(time_shift) - f(section_start))
            else
                pt = p0 + scale * (f(section_end) - f(section_start))
            end

            if isa(pd, AbstractDict)
                pd[pkey] = pt
            elseif isa(pd, AbstractVector)
                pd[Int(pkey)] = pt
            else
                # Try per-parameter setter on owning system if available
                if rss.owner !== nothing
                    try
                        DynamicalSystemsBase.set_parameter!(rss.owner, pkey, pt)
                    catch
                        # fall through to attempt setproperty!
                        try
                            setproperty!(pd, pkey, pt)
                        catch err
                            throw(ArgumentError("Cannot update parameter key $(pkey) for container type $(typeof(pd)). Provide a system with `set_parameter!` or use Dict/Vector parameter containers."))
                        end
                    end
                else
                    try
                        setproperty!(pd, pkey, pt)
                    catch err
                        throw(ArgumentError("Cannot update parameter key $(pkey) for container type $(typeof(pd)). Provide a RateSystem constructed from a system supporting `set_parameter!` or use Dict/Vector parameter containers."))
                    end
                end
            end
        end

        # If we updated the owning system per-key above, return its current parameters
        if rss.owner !== nothing
            try
                return DynamicalSystemsBase.current_parameters(rss.owner)
            catch
                # fallback: try to set all parameters at once
                try
                    DynamicalSystemsBase.set_parameters!(rss.owner, pd)
                    rss.pdummy = pd
                    return DynamicalSystemsBase.current_parameters(rss.owner)
                catch
                    rss.pdummy = pd
                    return pd
                end
            end
        end

        rss.pdummy = pd
        return pd
    end
end

# TODO: ensure all three key update functions follow the optional [, pidx] syntax
# for multi parameters.
"""
    set_forcing_scale!(rs::RateSystem, scale [, pidx])

Sets the amplitude (`forcing_scale`) of the forcing of the [`RateSystem`](@ref) to `scale`.
For multiple parameters, if `pidx` is given, change the forcing only corresponding to
the specified parameter, otherwise update the forcing scale of _all_ parameters to `scale.`
"""
function set_forcing_scale!(rs::RateSystem, scale; pidx=nothing)
    if pidx === nothing
        for k in keys(rs.forcing.forcers)
            rs.forcing.forcing_scale[k] = scale
        end
    else
        rs.forcing.forcing_scale[pidx] = scale
    end
    return rs
end

"""
$(TYPEDSIGNATURES)

Sets the duration (`forcing_duration`) of the forcing protocol applied to
the [`RateSystem`](@ref) `rs`.
"""
function set_forcing_duration!(rs::RateSystem, duration; pidx=nothing)
    if pidx === nothing
        for k in keys(rs.forcing.forcers)
            rs.forcing.forcing_duration[k] = duration
        end
    else
        rs.forcing.forcing_duration[pidx] = duration
    end
    return rs
end

"""
$(TYPEDSIGNATURES)

Sets the start time (`forcing_start_time`) of the forcing protocol applied to
the [`RateSystem`](@ref) `rs`.
"""
function set_forcing_start!(rs::RateSystem, start_time; pidx=nothing)
    if pidx === nothing
        for k in keys(rs.forcing.forcers)
            rs.forcing.forcing_start_time[k] = start_time
        end
    else
        rs.forcing.forcing_start_time[pidx] = start_time
    end
    return rs
end

"""
    current_parameters(rs::RateSystem, t)

Returns the parameter vector of a [`RateSystem`](@ref) at time `t` (in system time units).
Note: this function returns a copy of the original parameter container.
"""
function DynamicalSystemsBase.current_parameters(rs::RateSystem, t)
    p = deepcopy(current_parameters(rs.system))
    return p_modified(rs.forcing, p, t)
end

"""
$(TYPEDSIGNATURES)

Returns an autonomous dynamical system of type [`DynamicalSystemsBase.CoupledODEs`](@ref)
corresponding to the frozen system of the non-autonomous [`RateSystem`](@ref) `rs` at
time `t`.
"""
function frozen_system(rs::RateSystem, t)
    ds = CoupledODEs(rs.forcing.unforced_rule, current_state(rs), DynamicalSystemsBase.current_parameters(rs, t))
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
function SciMLBase.isinplace(rs::RateSystem)
    # Only advertise in-place if the underlying system is in-place and the
    # wrapped `unforced_rule` actually provides an in-place method for the
    # concrete argument types we will pass. This avoids promising in-place
    # semantics when we would need to mutate solver-owned parameter containers.
    try
        # If the wrapped system isn't in-place, we can't safely advertise it.
        if !SciMLBase.isinplace(rs.system)
            return false
        end

        # Obtain concrete argument types for the current system state/params.
        u = DynamicalSystemsBase.current_state(rs)
        p = DynamicalSystemsBase.current_parameters(rs)
        t = DynamicalSystemsBase.initial_time(rs)

        du_type = typeof(u)
        u_type = typeof(u)
        p_type = typeof(p)
        t_type = typeof(t)

        # Check whether an in-place method of the form f(du,u,p,t) exists
        return hasmethod(rs.forcing.unforced_rule, Tuple{du_type, u_type, p_type, t_type})
    catch
        return false
    end
end
StateSpaceSets.dimension(rs::RateSystem) = StateSpaceSets.dimension(rs.system)
DynamicalSystemsBase.isdeterministic(rs::RateSystem) = true # TODO: only true if CoupledODEs

(rs::RateSystem)(t::Real) = rs.system(t)
