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

Mutable container storing the information required to construct and operate a
`RateSystem` implementation.

Fields include the underlying autonomous `unforced_rule`, a mapping of
parameter keys to `ForcingProfile` instances (`forcers`), per-key timing and
scale maps (`forcing_start_time`, `forcing_duration`, `forcing_scale`), the
initial autonomous parameter values (`p0`), an initial time `t0`, a cached
dummy parameter container `pdummy` (used when the concrete parameter container
is not a `Dict` or `Vector`), and a reference to the owning system `owner`.

Calling semantics:
- `(::RateSystemSpecs)(u, p, t)` — out-of-place call: returns the computed
    derivative `du` (the nonautonomous drift) using a modified parameter container appropriate for time
    `t`, where `p` is the parameter container of the underlying autonomous system.
- `(::RateSystemSpecs)(du, u, p, t)` — in-place call: attempts to call an
    in-place `unforced_rule(du, u, p, t)` and, if not available, falls back to
    the out-of-place call and copies the result into `du`.

The `RateSystemSpecs` ensures parameter values are adjusted according to the
configured forcing profiles before delegating to the underlying autonomous rule.
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
"""
    initial_parameter_map(rss::RateSystemSpecs)

Return the internal mapping (`Dict`) of parameter keys to the initial (autonomous)
parameter values stored in `rss.p0`.
"""
function initial_parameter_map(rss::RateSystemSpecs)
    return getfield(rss, :p0)
end

"""
    initial_parameter(rss::RateSystemSpecs, pkey)

Return the initial autonomous parameter value associated with `pkey`.
"""
function initial_parameter(rss::RateSystemSpecs, pkey)
    return getfield(rss, :p0)[pkey]
end

"""
    initial_parameter_value(rss::RateSystemSpecs)

Return the single initial autonomous parameter value when exactly one forcer
is present. Throws `ArgumentError` if multiple forcers are configured.
"""
function initial_parameter_value(rss::RateSystemSpecs)
    p0map = getfield(rss, :p0)
    if length(getfield(rss, :forcers)) == 1
        return first(values(p0map))
    else
        throw(ArgumentError("initial_parameter_value(rss) only valid when exactly one forcer is present"))
    end
end

# RateSystem-level wrappers:
"""
    parameters(sys, args...; kw...)

Convenience wrapper delegating to `current_parameters(sys, ...)` for API
compatibility with earlier code and callers expecting a `parameters` helper.
"""
function parameters(sys, args...; kw...)
    return current_parameters(sys, args...; kw...)
end

"""
        RateSystem <: ContinuousTimeDynamicalSystem
        RateSystem(ds::ContinuousTimeDynamicalSystem, forcing_profiles::Dict; kw...)

Construct a non-autonomous dynamical system by applying time-dependent parametric
forcing protocols to an underlying autonomous continuous-time dynamical system `ds`.

`forcing_profiles` must be a `Dict` mapping parameter indices (anything valid for
`set_parameter!`) to `ForcingProfile` instances; each entry defines the functional 
form of how the corresponding parameter evolves over time.

Keyword arguments
- `forcing_start_time` (default: `nothing`) — if `nothing`, each forcer's start
    time is set to `t0` (the system initial time). You may supply an `AbstractDict`
    mapping keys to start times, or a single scalar value which will be applied to 
    all forcers, giving the time(s) for which the ramping of the corresponding 
    parameters starts.
- `forcing_duration` (default: `nothing`) — if `nothing`, each forcer's
    duration defaults to `fp.interval[2] - fp.interval[1]` for that profile. Can
    be an `AbstractDict` or a scalar applied to all forcers, giving the duration of 
    the parameter ramping (in system time units).
- `forcing_scale` (default: `nothing`) — if `nothing`, defaults to `1.0` for
    each forcer. Can be an `AbstractDict` or a scalar applied to all forcers. 
    Acts as amplitude multiplication factor of the forcing protocol(s).
- `t0` (default: `initial_time(ds)`) — initial time of the `RateSystem`.

Description
The profile `interval` defines the domain of the forcing function; when applied
to the system the profile is rescaled in system time units using the
configured `start` and `duration` values - this allow changing the rate of the 
parameter ramping. Before a forcer's `start` time the parameter equals its initial 
autonomous value; during the forcing interval it follows the rescaled profile (multiplied 
by the corresponding `forcing_scale` factor); after the interval the parameter is frozen 
at its final value.

Multiple parameters
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

    # Wrap the modified (non-autonomous) drift in an appropriate Coupled* wrapper.
    # If the input system is stochastic, construct a CoupledSDEs wrapper and
    # attempt to preserve common SDE-related keyword settings (e.g. `g`,
    # `noise_prototype`, `noise_strength`, `covariance`, `diffeq`). Otherwise
    # fall back to a CoupledODEs wrapper.
    if ds isa CoupledSDEs
        kw = (;)
        if hasproperty(ds, :g)
            kw = merge(kw, (; g = getproperty(ds, :g)))
        end
        if hasproperty(ds, :noise_prototype)
            kw = merge(kw, (; noise_prototype = getproperty(ds, :noise_prototype)))
        end
        if hasproperty(ds, :noise_strength)
            kw = merge(kw, (; noise_strength = getproperty(ds, :noise_strength)))
        end
        if hasproperty(ds, :covariance)
            kw = merge(kw, (; covariance = getproperty(ds, :covariance)))
        end
        if hasproperty(ds, :diffeq)
            kw = merge(kw, (; diffeq = getproperty(ds, :diffeq)))
        end

        # preserve initial time in wrapper
        kw = merge(kw, (; t0 = t0 ))

        newds = CoupledSDEs(rss, current_state(ds), deepcopy(current_parameters(ds)); kw...)
    else
        newds = CoupledODEs(rss, current_state(ds), deepcopy(current_parameters(ds)); t0)
    end

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

"""
        p_modified(rss::RateSystemSpecs, p, t)

Compute and return a parameter container appropriate for time `t` by applying
the forcing profiles configured in `rss.forcers` to the provided parameter
container `p`.

Behavior:
- If `p` is an `AbstractDict` or `AbstractVector`, it is updated in-place and
    returned (this avoids allocations when possible).
- For arbitrary parameter container types, a copy of `rss.pdummy` is created and
    updated. When `rss.owner` is available, `DynamicalSystemsBase.set_parameter!`
    (per-key) or `set_parameters!` (bulk) is used where supported; otherwise
    `setproperty!` is attempted on the copy.

The returned container has the same general shape/type expected by the
underlying system.
"""
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
