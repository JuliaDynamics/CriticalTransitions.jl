"""
    RateSystemSpecs <: Function

A mutable data structure storing information needed to construct and modify a
[`RateSystem`](@ref).

## Fields
$(TYPEDFIELDS)

Call signature: `(::RateSystemSpecs)(u, p, t)` for out-of-place and
`(::RateSystemSpecs)(du, u, p, t)` for in-place dynamical systems.
"""
mutable struct RateSystemSpecs{S,K,T,P} <: Function
    "Underlying autonomous system"
    unforced_system::S
    "Mapping parameter index => ForcingProfile"
    forcing_profile::AbstractDict{K,ForcingProfile}
    "Mapping parameter index => forcing start time"
    forcing_start_time::AbstractDict{K,T}
    "Mapping parameter index => forcing duration"
    forcing_duration::AbstractDict{K,T}
    "Mapping parameter index => forcing scale"
    forcing_scale::AbstractDict{K,T}
    "Placeholder parameter container"
    pdummy::P
    "Initial time (of system initiation)"
    t0::T
end

"""
    RateSystem <: ContinuousTimeDynamicalSystem
    RateSystem(ds::ContinuousTimeDynamicalSystem, forcing_profile::AbstractDict; kw...)

Construct a non-autonomous dynamical system by applying time-dependent parametric
forcing protocols to an underlying autonomous continuous-time dynamical system `ds`.

`forcing_profile` must be a `Dict` mapping parameter indices (anything valid for
`set_parameter!`) to `ForcingProfile` instances; each entry defines the functional 
form of how the corresponding parameter evolves over time.

## Keyword arguments
- `forcing_start_time` (default: `initial_time(ds)`) — start time(s) for the parameter
    ramping. You may supply an `AbstractDict` mapping keys to start times, or a scalar
    value which will be applied to all forcing profiles.
- `forcing_duration` (default: `1.0`) — duration of the parameter ramping (in system
    time units). Can be an `AbstractDict` mapping keys to durations or a scalar
    applied to all forcing profiles.
- `forcing_scale` (default: `1.0`) — amplitude multiplicative factor of the forcing
    profile. Can be an `AbstractDict` mapping keys to scales or a scalar applied to
    all forcing profiles.
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
    specs::RateSystemSpecs{R}
end

function RateSystem(
    ds::ContinuousTimeDynamicalSystem,
    forcing_profile::AbstractDict;
    forcing_start_time=initial_time(ds),
    forcing_duration=1.0,
    forcing_scale=1.0,
    t0=initial_time(ds),
    )

    if isempty(forcing_profile)
        throw(ArgumentError("`forcing_profile` must be a dictionary containing at least one entry"))
    end

    forcing_kw = Dict("start_time" => forcing_start_time, "duration" => forcing_duration,
        "scale" => forcing_scale)
    
    for j in keys(forcing_kw)
        if forcing_kw[j] isa Real
            forcing_kw[j] = Dict(k => forcing_kw[j] for (k, _) in forcing_profile)
        else
            @assert keys(forcing_kw[j]) == keys(forcing_profile) |
            "Dictionary keys of forcing_$(j) and forcing_profile must be identical."
        end
    end

    p0 = deepcopy(current_parameters(ds))

    rss = RateSystemSpecs{R,K,T,P}(
        ds,
        forcing_profile,
        forcing_kw["start_time"],
        forcing_kw["duration"],
        forcing_kw["scale"],
        p0,
        t0
    )

    # preserve CoupledSDEs properties
    if ds isa CoupledSDEs
        IIP = SciMLBase.isinplace(ds)
        prob = referrenced_sciml_prob(ds)
        new_prob = SDEProblem{IIP}(
            rss, prob.g, current_state(ds), (t0, prob.tspan[2]), p0;
            noise_rate_prototype = prob.noise_rate_prototype,
            noise = prob.noise,
        )
        system = CoupledSDEs(new_prob, ds.diffeq, ds.noise_type)
    elseif ds isa CoupledODEs
        system = CoupledODEs(rss, current_state(ds), p0; t0)
    else
        error("A RateSystem can only be constructed from a CoupledODEs or CoupledSDEs.")
    end

    return RateSystem(system, rss)
end

RateSystem(ds::ContinuousTimeDynamicalSystem, forcing_profile::ForcingProfile,
    pkey; kw...) = RateSystem(ds, Dict(pkey => forcing_profile); kw...)

# Out-of-place
function (rss::RateSystemSpecs)(u, p, t)
    update_parameters!(rss, t)
    return dynamic_rule(rss.unforced_system)(u, rss.pdummy, rss.t0)
end

# In-place
function (rss::RateSystemSpecs)(du, u, p, t)
    update_parameters!(rss, t)
    return dynamic_rule(rss.unforced_system)(du, u, rss.pdummy, rss.t0)
end

function update_parameters!(rss::RateSystemSpecs, t::Real)
    p_dummy = rss.pdummy
    ds_dummy = rss.unforced_system

    for (pkey, profile) in rss.forcing_profile
        p_old = current_parameter(ds, pkey)
        f = profile.profile

        section_start = profile.interval[1]
        section_end = profile.interval[2]
        # Make the function piecewise constant
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
    set_parameters!(ds_dummy, p_dummy, rss.pdummy)
    return nothing
end

"""
    parameters(rs::RateSystem, t)

Returns the parameter container of a [`RateSystem`](@ref) at time `t` (in system time
units).
"""
function parameters(rs::RateSystem, t)
    update_parameters!(rs.specs, t)
    return rs.specs.pdummy
end

"""
    parameter(rs::RateSystem, t, pkey)

Returns the parameter with key `pkey` of a [`RateSystem`](@ref) at time `t` (in system
time units).
"""
function parameter(rs::RateSystem, t, pkey)
    update_parameters!(rs.specs, t)
    return rs.specs.pdummy[pkey]
end

"""
    set_forcing_scale!(rs::RateSystem, scale; pkey=nothing)

Sets the amplitude (`forcing_scale`) of the forcing of the [`RateSystem`](@ref) to `scale`.
If the keyword `pkey` is provided, only the corresponding parameter's forcing scale is
changed; otherwise the forcing scale of all parameters is updated to `scale`.
"""
function set_forcing_scale!(rs::RateSystem, scale; pkey=nothing)
    if isnothing(pkey)
        for k in keys(rs.specs.forcing_profile)
            rs.specs.forcing_scale[k] = scale
        end
    else
        rs.specs.forcing_scale[pkey] = scale
    end
    return rs
end

"""
$(TYPEDSIGNATURES)

Sets the duration (`forcing_duration`) of the forcing protocol applied to
the [`RateSystem`](@ref) `rs`. If the optional keyword `pkey` is provided, only the
duration for that parameter is changed; otherwise all parameters are updated.
"""
function set_forcing_duration!(rs::RateSystem, duration; pkey=nothing)
    if isnothing(pkey)
        for k in keys(rs.specs.forcing_profile)
            rs.specs.forcing_duration[k] = duration
        end
    else
        rs.specs.forcing_duration[pkey] = duration
    end
    return rs
end

"""
$(TYPEDSIGNATURES)

Sets the start time (`forcing_start_time`) of the forcing protocol applied to
the [`RateSystem`](@ref) `rs`. If the optional keyword `pkey` is provided, only the
start time for that parameter is changed; otherwise all parameters are updated.
"""
function set_forcing_start!(rs::RateSystem, start_time; pkey=nothing)
    if isnothing(pkey)
        for k in keys(rs.specs.forcing_profile)
            rs.specs.forcing_start_time[k] = start_time
        end
    else
        rs.specs.forcing_start_time[pkey] = start_time
    end
    return rs
end

"""
$(TYPEDSIGNATURES)

Returns an autonomous dynamical system corresponding to the frozen system of the
non-autonomous [`RateSystem`](@ref) `rs` at time `t`.
"""
function frozen_system(rs::RateSystem, t)
    p = parameters(rs, t)
    # preserve CoupledSDEs properties
    if rs isa CoupledSDEs
        IIP = SciMLBase.isinplace(rs)
        prob = referrenced_sciml_prob(rs)
        new_prob = SDEProblem{IIP}(
            dynamic_rule(rs.specs.unforced_system), prob.g, current_state(rs),
            (t0, prob.tspan[2]), p;
            noise_rate_prototype = prob.noise_rate_prototype,
            noise = prob.noise,
        )
        return CoupledSDEs(new_prob, rs.diffeq, rs.noise_type)
    else
        return CoupledODEs(dynamic_rule(rs.specs.unforced_system), current_state(rs), p)
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
SciMLBase.isinplace(rs::RateSystem) = SciMLBase.isinplace(rs.system)
DynamicalSystemsBase.set_state!(rs::RateSystem, u::AbstractArray) = (DynamicalSystemsBase.set_state!(rs.system, u); rs)

# SDE-specific delegations
noise_process(rs::RateSystem) = noise_process(rs.system)

solver(rs::RateSystem) = solver(rs.system)

DynamicalSystemsBase.covariance_matrix(rs::RateSystem, args...; kw...) = covariance_matrix(rs.system, args...; kw...)

DynamicalSystemsBase.diffusion_matrix(rs::RateSystem, args...; kw...) = diffusion_matrix(rs.system, args...; kw...)

StateSpaceSets.dimension(rs::RateSystem) = StateSpaceSets.dimension(rs.system)
DynamicalSystemsBase.isdeterministic(rs::RateSystem) = isa(rs.system, CoupledODEs)
(rs::RateSystem)(t::Real) = rs.system(t)
