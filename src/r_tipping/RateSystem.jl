"""
    RateSystemSpecs <: Function

A mutable data structure storing information needed to construct and modify a
[`RateSystem`](@ref).

## Fields
$(FIELDS)

Call signature: `(::RateSystemSpecs)(u, p, t)` for out-of-place and
`(::RateSystemSpecs)(du, u, p, t)` for in-place dynamical systems.
"""
mutable struct RateSystemSpecs{R,K,T,P} <: Function
    "Dynamic rule of the underlying autonomous system"
    unforced_rule::R
    "Mapping parameter index => ForcingProfile"
    forcing_profile::Dict{K,ForcingProfile}
    "Mapping parameter index => forcing start time"
    forcing_start_time::Dict{K,T}
    "Mapping parameter index => forcing duration"
    forcing_duration::Dict{K,T}
    "Mapping parameter index => forcing scale"
    forcing_scale::Dict{K,T}
    "Parameter container of the autonomous system (at t0)"
    p0::P
    "Placeholder parameter container (used for update_parameters!)"
    pdummy::P
    "Initial time (of system initiation)"
    t0::T
end

"""
    RateSystem <: ContinuousTimeDynamicalSystem
    RateSystem(ds::ContinuousTimeDynamicalSystem, forcing_profile::Dict; kw...)

Construct a non-autonomous dynamical system by applying time-dependent parametric
forcing protocols to an underlying autonomous continuous-time dynamical system `ds`.

`forcing_profile` must be a `Dict` mapping parameter indices (anything valid for
`set_parameter!`) to `ForcingProfile` instances; each entry defines the functional 
form of how the corresponding parameter evolves over time.

## Keyword arguments
- `forcing_start_time` (default: `nothing`) — if `nothing`, each forcing_profile's start
    time is set to `t0` (the system initial time). You may supply an `AbstractDict`
    mapping keys to start times, or a single scalar value which will be applied to 
    all forcing_profile, giving the time(s) for which the ramping of the corresponding 
    parameters starts.
- `forcing_duration` (default: `nothing`) — if `nothing`, each forcing_profile's
    duration defaults to `fp.interval[2] - fp.interval[1]` for that profile. Can
    be an `AbstractDict` or a scalar applied to all forcing_profile, giving the duration of 
    the parameter ramping (in system time units).
- `forcing_scale` (default: `nothing`) — if `nothing`, defaults to `1.0` for
    each forcing_profile. Can be an `AbstractDict` or a scalar applied to all forcing_profile. 
    Acts as amplitude multiplication factor of the forcing protocol(s).
- `t0` (default: `initial_time(ds)`) — initial time of the `RateSystem`.

## Description
The profile `interval` defines the domain of the forcing function; when applied
to the system the profile is rescaled in system time units using the
configured `start` and `duration` values - this allow changing the rate of the 
parameter ramping. Before a forcing_profile's `start` time the parameter equals its initial 
autonomous value; during the forcing interval it follows the rescaled profile (multiplied 
by the corresponding `forcing_scale` factor); after the interval the parameter is frozen 
at its final value.

## Multiple parameters
Pass a `Dict` with one `ForcingProfile` per parameter to force multiple
parameters simultaneously. The `forcing_start_time`, `forcing_duration`, and
`forcing_scale` keywords accept the same flexibility (per-key `Dict` or scalar
applied to all keys) in this mode.
"""
struct RateSystem{S,R} <: ContinuousTimeDynamicalSystem
    "Non-autonomous continuous-time dynamical system"
    system::S
    "Data structure representing the system and forcing specifications in system units"
    forcing::RateSystemSpecs{R}
end

function RateSystem(
    ds::ContinuousTimeDynamicalSystem,
    forcing_profile::Dict;
    forcing_start_time::Dict = Dict(),
    forcing_duration::Dict = Dict(),
    forcing_scale::Dict = Dict(),
    t0=initial_time(ds),
    )

    if isempty(forcing_profile)
        throw(ArgumentError("`forcing_profile` dictionary must contain at least one entry"))
    end

    if isempty(forcing_start_time)
        for (k, _) in forcing_profile
            forcing_start_time[k] = initial_time(ds)
        end
    end

    if isempty(forcing_duration)
        for (k, profile) in forcing_profile
            forcing_duration[k] = (profile.interval[2] - profile.interval[1])
        end
    end

    if isempty(forcing_scale)
        for (k, _) in forcing_profile
            forcing_scale[k] = 1.0
        end
    end

    p0 = deepcopy(current_parameters(ds))

    rss = RateSystemSpecs{R,K,T,P}(
        dynamic_rule(ds),
        forcing_profile,
        forcing_start_time,
        forcing_duration,
        forcing_scale,
        p0,
        p0,
        t0
    )

    # preserve CoupledSDEs properties
    if ds isa CoupledSDEs
        kw = (;)
        hasproperty(ds, :g) ? kw = merge(kw, (; g = getproperty(ds, :g))) : nothing
        hasproperty(ds, :noise_prototype) ? kw = merge(kw,
            (; noise_prototype = getproperty(ds, :noise_prototype))) : nothing
        hasproperty(ds, :noise_strength) ? kw = merge(kw,
            (; noise_strength = getproperty(ds, :noise_strength))) : nothing
        hasproperty(ds, :covariance) ? kw = merge(kw,
            (; covariance = getproperty(ds, :covariance))) : nothing
        hasproperty(ds, :diffeq) ? kw = merge(kw,
            (; diffeq = getproperty(ds, :diffeq))) : nothing

        # preserve initial time in wrapper
        kw = merge(kw, (; t0 = t0 ))

        system = CoupledSDEs(rss, current_state(ds), p0; kw...)
    elseif ds isa CoupledODEs
        system = CoupledODEs(rss, current_state(ds), p0; t0)
    else
        error("A RateSystem can only be constructed from a CoupledODEs or CoupledSDEs.")
    end

    return RateSystem(system, rss)
end

function RateSystem(ds::ContinuousTimeDynamicalSystem, forcing_profile::ForcingProfile,
    pkey;
    forcing_start_time::Real = 0.0, 
    forcing_duration::Real = 1.0,
    forcing_scale::Real = 1.0,
    kwargs...)

    return RateSystem(ds, Dict(pkey => forcing_profile);
        forcing_start_time = Dict(pkey => forcing_start_time),
        forcing_duration = Dict(pkey => forcing_duration),
        forcing_scale = Dict(pkey => forcing_scale),
        kwargs...)
end

# Out-of-place
function (rss::RateSystemSpecs)(u, p, t)
    update_parameters!(rss, t)
    return rss.unforced_rule(u, rss.pdummy, rss.t0)
end

# In-place
function (rss::RateSystemSpecs)(du, u, p, t)
    update_parameters!(rss, t)
    return rss.unforced_rule(du, u, rss.pdummy, rss.t0)
end

function update_parameters!(rss::RateSystemSpecs, t::Real)
    pdummy = rss.pdummy

    for (pkey, profile) in rss.forcing_profile
        p_old = current_parameter(ds, pkey)
        f = profile.profile

        section_start = profile.interval[1]
        section_end = profile.interval[2]
        # Make the function piecewise constant with range [p0, p0+forcing_scale]
        if t > rss.forcing_start_time[pkey]
            if t < rss.forcing_start_time[pkey] + rss.forcing_duration[pkey]
                # Stretch/squeeze forcing to the correct time units
                time_shift =
                    ((section_end - section_start) / rss.forcing_duration[pkey]) *
                    (t - rss.forcing_start_time[pkey]) + section_start
                p_new = p_old + rss.forcing_scale[pkey] * 
                    (f(time_shift) - f(section_start))
            else
                p_new = p_old + rss.forcing_scale[pkey] * 
                    (f(section_end) - f(section_start))
            end
        else
            p_new = p_old
        end
        pdummy[pkey] = p_new
    end
    set_parameters!(rss, pdummy, rss.pdummy)
    return nothing
end

"""
    parameters(rs::RateSystem, t)

Returns the parameter container of a [`RateSystem`](@ref) at time `t` (in system time
units).
"""
function parameters(rs::RateSystem, t)
    update_parameters!(rs.forcing, t)
    return rs.forcing.pdummy
end

"""
    parameter(rs::RateSystem, t, pkey)

Returns the parameter with key `pkey` of a [`RateSystem`](@ref) at time `t` (in system
time units).
"""
function parameter(rs::RateSystem, t, pkey)
    update_parameters!(rs.forcing, t)
    return rs.forcing.pdummy[pkey]
end

"""
    set_forcing_scale!(rs::RateSystem, scale [, pkey])

Sets the amplitude (`forcing_scale`) of the forcing of the [`RateSystem`](@ref) to `scale`.
For multiple parameters, if `pkey` is given, change the forcing only corresponding to
the specified parameter, otherwise update the forcing scale of _all_ parameters to `scale.`
"""
function set_forcing_scale!(rs::RateSystem, scale; pkey=nothing)
    if isnothing(pkey)
        for k in keys(rs.forcing.forcing_profile)
            rs.forcing.forcing_scale[k] = scale
        end
    else
        rs.forcing.forcing_scale[pkey] = scale
    end
    return rs
end

"""
$(TYPEDSIGNATURES)

Sets the duration (`forcing_duration`) of the forcing protocol applied to
the [`RateSystem`](@ref) `rs`.
"""
function set_forcing_duration!(rs::RateSystem, duration; pkey=nothing)
    if isnothing(pkey)
        for k in keys(rs.forcing.forcing_profile)
            rs.forcing.forcing_duration[k] = duration
        end
    else
        rs.forcing.forcing_duration[pkey] = duration
    end
    return rs
end

"""
$(TYPEDSIGNATURES)

Sets the start time (`forcing_start_time`) of the forcing protocol applied to
the [`RateSystem`](@ref) `rs`.
"""
function set_forcing_start!(rs::RateSystem, start_time; pkey=nothing)
    if isnothing(pkey)
        for k in keys(rs.forcing.forcing_profile)
            rs.forcing.forcing_start_time[k] = start_time
        end
    else
        rs.forcing.forcing_start_time[pkey] = start_time
    end
    return rs
end

"""
$(TYPEDSIGNATURES)

Returns an autonomous dynamical system corresponding to the frozen system of the
non-autonomous [`RateSystem`](@ref) `rs` at time `t`. If the wrapped system is an
SDE, a `CoupledSDEs` is returned (attempting to preserve the original SDE's
noise/integrator settings); otherwise a `CoupledODEs` is returned.
"""
function frozen_system(rs::RateSystem, t)
    p = DynamicalSystemsBase.current_parameters(rs, t)
    if rs.system isa CoupledSDEs
        kw = (;)
        if hasproperty(rs.system, :g)
            kw = merge(kw, (; g = getproperty(rs.system, :g)))
        end
        if hasproperty(rs.system, :noise_prototype)
            kw = merge(kw, (; noise_prototype = getproperty(rs.system, :noise_prototype)))
        end
        if hasproperty(rs.system, :noise_strength)
            kw = merge(kw, (; noise_strength = getproperty(rs.system, :noise_strength)))
        end
        if hasproperty(rs.system, :covariance)
            kw = merge(kw, (; covariance = getproperty(rs.system, :covariance)))
        end
        if hasproperty(rs.system, :diffeq)
            kw = merge(kw, (; diffeq = getproperty(rs.system, :diffeq)))
        end

        # preserve the RateSystem initial time for the frozen wrapper
        kw = merge(kw, (; t0 = rs.forcing.t0 ))

        return CoupledSDEs(rs.forcing.unforced_rule, current_state(rs), p; kw...)
    else
        return CoupledODEs(rs.forcing.unforced_rule, current_state(rs), p)
    end
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
DynamicalSystemsBase.set_state!(rs::RateSystem, u::AbstractArray) = (DynamicalSystemsBase.set_state!(rs.system, u); rs)

# SDE-specific delegations: forward stochastic helpers to the wrapped system
function noise_process(rs::RateSystem)
    return noise_process(rs.system)
end

function solver(rs::RateSystem)
    return solver(rs.system)
end

function StochasticSystemsBase.covariance_matrix(rs::RateSystem, args...; kw...)
    return StochasticSystemsBase.covariance_matrix(rs.system, args...; kw...)
end

function StochasticSystemsBase.diffusion_matrix(rs::RateSystem, args...; kw...)
    return StochasticSystemsBase.diffusion_matrix(rs.system, args...; kw...)
end

StateSpaceSets.dimension(rs::RateSystem) = StateSpaceSets.dimension(rs.system)
DynamicalSystemsBase.isdeterministic(rs::RateSystem) = isa(rs.system, CoupledODEs)
(rs::RateSystem)(t::Real) = rs.system(t)
