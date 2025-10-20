mutable struct NonautonomousDynamicRule{R,F,P,T} <: Function
    unforced_rule::R
    pidx::P
    fp::ForcingProfile{F,T}
    forcing_start_time::T
    forcing_length::T
    forcing_scale::T
end

"""
    RateSystem <: ContinuousTimeDynamicalSystem

A non-autonomous dynamical system obtained by applying a time-dependent parametric
forcing protocol `profile` to an underlying autonomous continuous-time dynamical system
`ds`.

The time-dependent forcing protocol `profile` describes the evolution of the parameter defined
by `pidx` over a time interval `interval`. The parameter `p` that is swept is defined by
`current_parameters(ds)[pidx]`.

The forcing protocol can be specified directly via a `[ForcingProfile](@def)` referred to as
`fp = ForcingProfile(profile, interval)`. The forcing protocol can be adjusted by using the
keyword arguments `forcing_start_time`, `forcing_length`, and `forcing_scale` when
constructing the `RateSystem`. This allows to repurpose the same forcing protocol `fp`.
"""
struct RateSystem{S,F,P,T} <: ContinuousTimeDynamicalSystem
    "Underlying continuous-time dynamical system"
    system::S
    "The index of the parameter being forced"
    forcing::NonautonomousDynamicRule{F,P,T}
end

"""
    RateSystem(ds, fp::ForcingProfile, pidx; kwargs...)

## Keyword arguments
- `forcing_start_time = interval[1]`: Time when parameter shift starts
  (before this, the resulting system will be autonomous)
- `forcing_length = interval[2] - interval[1]`: Time-interval over which `profile(interval)` is
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
    forcing_length::T=(fp.interval[2] - fp.interval[1]),
    forcing_scale::T=1.0,
    t0::T=initial_time(ds),
) where {P,T<:Real}
    (forcing_start_time >= t0) ||
        throw("The forcing cannot start before the system's initial time `t0`
              but your forcing_start_time ($(forcing_start_time)) < t0 ($(t0)).")
    _ds = deepcopy(ds)
    forcing_config = NonautonomousDynamicRule(
        dynamic_rule(_ds), pidx, fp, forcing_start_time, forcing_length, forcing_scale
    )
    system = CoupledODEs(forcing_config, current_state(_ds), current_parameters(_ds); t0)
    return RateSystem(system, forcing_config)
end

function (ndr::NonautonomousDynamicRule)(u, p, t)
    # this should be rewritten with the Callback mechanism from SciML
    p_i = p[ndr.pidx]
    p_at_t = p_modified(t, ndr, p_i)
    if p isa Union{AbstractArray,AbstractDict}
        setindex!(p, p_at_t, ndr.pidx)
    else
        setproperty!(p, ndr.pidx, p_at_t)
    end
    return ndr.unforced_rule(u, p, t)
end

"""
Creates a piecewise constant function in alignment with the entries of the ForcingProfile
and the parameter value of the underlying autonomous system
"""
function p_modified(t::Real, ndr::NonautonomousDynamicRule, p0)
    f = ndr.fp.profile
    section_start = ndr.fp.interval[1]
    section_end = ndr.fp.interval[2]
    # making the function piecewise constant with range [p0,p0+forcing_scale]
    q = if t ≤ ndr.forcing_start_time
        return p0
    else
        if t < ndr.forcing_start_time + ndr.forcing_length
            # performing the time shift corresponding to stretching/squeezing
            time_shift =
                ((section_end - section_start) / ndr.forcing_length) *
                (t - ndr.forcing_start_time) + section_start
            return p0 + ndr.forcing_scale*(f(time_shift) - f(section_start))
        else
            return p0 + ndr.forcing_scale*(f(section_end) - f(section_start))
        end
    end
    return q
end

function set_forcing_scale!(rs::RateSystem, scale)
    rs.ndr.forcing_scale = scale
    return rs
end

function set_forcing_length!(rs::RateSystem, length)
    rs.ndr.forcing_length = length
    return rs
end

function set_forcing_start!(rs::RateSystem, start)
    rs.ndr.forcing_start = start
    return rs
end

function get_forcing(rs, t)
    p = initial_parameters(rs)
    p[rs.ndr.pidx] = p_modified(t, rs.ndr, p)
    return p
end

function frozen_system(rs::RateSystem, t)
    ds = CoupledODEs(rs.ndr.unforced_rule, current_state(ds), get_forcing(rs, t))
    return ds
end

# Extensions
for f in (
    :initial_state,
    :initial_parameters,
    :current_parameters,
    :dynamic_rule,
    :current_time,
    :current_state,
    :initial_time,
    :successful_step,
    :set_parameter!,
    :trajectory,
) # all api functions here
    @eval DynamicalSystemsBase.$(f)(rs::RateSystem, args...; kw...) = $(f)(
        rs.system, args...; kw...
    )
end

function DynamicalSystemsBase.set_state!(rs::RateSystem, u::AbstractArray)
    return (DynamicalSystemsBase.set_state!(rs.system, u); rs)
end
SciMLBase.step!(rs::RateSystem, args...) = (SciMLBase.step!(rs.system, args...); rs)
SciMLBase.isinplace(rs::RateSystem) = SciMLBase.isinplace(rs.system)
StateSpaceSets.dimension(rs::RateSystem) = StateSpaceSets.dimension(rs.system)
DynamicalSystemsBase.isdeterministic(rs::RateSystem) = true # TODO: only true if CoupledODEs

(rs::RateSystem)(t::Real) = rs.system(t)
