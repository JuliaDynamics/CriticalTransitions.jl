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

A non-autonomous dynamical system constructed by applying a time-dependent parametric
forcing protocol ([`ForcingProfile`](@ref)) to an underlying autonomous continuous-time
dynamical system `ds`.
"""
struct RateSystem{S,F,P,T} <: ContinuousTimeDynamicalSystem
    "Non-autonomous continuous-time dynamical system"
    system::S
    "Data structure representing the system and forcing specifications in system units"
    forcing::RateSystemSpecs{F,P,T}
end

"""
    RateSystem(ds::ContinuousTimeDynamicalSystem, fp::ForcingProfile, pidx; kwargs...)

Constructs a non-autonomous dynamical system from a given autonomous 
continous-time dynamical system `ds` by applying a time-dependent parametric forcing 
via a [`ForcingProfile`](@ref).

The `ForcingProfile`, applied to the parameter with index `pidx`, describes the functional
form of the parameter evolution over a given time interval. This time interval is defined
in system time units by its start time (`forcing_start_time`) and duration
(`forcing_duration`). The forcing profile can be scaled in magnitude by the factor
`forcing_scale`.

Before `t = forcing_start_time`, the system parameters correspond to that of the
underlying autonomous system `ds`. At the end of the forcing window
(`t = forcing_start_time + forcing_duration`), the parameters are frozen at their current
value and thereafter the system is again autonomous. 

## Keyword arguments
- `forcing_start_time = interval[1]`: Time when parameter shift starts
  (before this, the resulting system will be autonomous)
- `forcing_duration = interval[2] - interval[1]`: Time-interval over which `profile(interval)` is
  spread out (for window_length > interval[2] - interval[1]) or squeezed into
  (for window_length < interval[2] - interval[1])
- `forcing_scale = 1.0`: Amplitude of the protocol.
- `t0 = 0.0`: Initial time of the resulting non-autonomous system
"""
function RateSystem(
    ds::ContinuousTimeDynamicalSystem,
    fp::ForcingProfile,
    pidx::P;
    forcing_start_time::T=initial_time(ds),
    forcing_duration::T=(fp.interval[2] - fp.interval[1]),
    forcing_scale::T=1.0,
    t0::T=initial_time(ds),
) where {P,T<:Real}
    (forcing_start_time >= t0) ||
        throw("The forcing cannot start before the system's initial time `t0`
              but your forcing_start_time ($(forcing_start_time)) < t0 ($(t0)).")

    # TODO: Make p0 a vector
    p0 = current_parameter(ds, pidx)
    rss = RateSystemSpecs(
        dynamic_rule(ds),
        pidx,
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
    if p isa Union{AbstractArray,AbstractDict}
        setindex!(p, p_at_t, rss.pidx)
    else
        setproperty!(p, rss.pidx, p_at_t)
    end
    # TODO: This doesn't work if the dynamical system is in-place...
    # We need to multi-dispatch one more version of `(rss::)` function with 4 arguments
    return rss.unforced_rule(u, p, rss.t0)
end

"""
Creates a piecewise constant function in alignment with the entries of the ForcingProfile
and the parameter value of the underlying autonomous system
"""
function p_modified(rss::RateSystemSpecs, t::Real)
    # TODO: Make `p0` vector, and make the function re-write a dummy vector of time dependent parameters
    p0 = rss.p0
    f = rss.fp.profile
    section_start = rss.fp.interval[1]
    section_end = rss.fp.interval[2]
    # making the function piecewise constant with range [p0,p0+forcing_scale]
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

function set_forcing_scale!(rs::RateSystem, scale)
    rs.forcing.forcing_scale = scale
    return rs
end

function set_forcing_duration!(rs::RateSystem, length)
    rs.forcing.forcing_duration = length
    return rs
end

function set_forcing_start!(rs::RateSystem, start)
    rs.forcing.forcing_start = start
    return rs
end

function get_forcing(rs, t)
    p = deepcopy(current_parameters(rs))
    # TODO: Doesn't work for vector parameters
    p[rs.forcing.pidx] = p_modified(rs.forcing, t)
    return p
end

function frozen_system(rs::RateSystem, t)
    ds = CoupledODEs(rs.forcing.unforced_rule, current_state(ds), get_forcing(rs, t))
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
