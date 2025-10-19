"""
    RateSystem <: ContinuousTimeDynamicalSystem
    RateSystem(ds, pfunc::Function, ptspan::Tuple, pidx; kwargs...)
    RateSystem(ds, rf::RateFunction, pidx; kwargs...)

A non-autonomous dynamical system obtained by applying a time-dependent parametric
forcing protocol `pfunc` to an underlying autonomous continuous-time dynamical system
`ds <: ContinuousTimeDynamicalSystem`.

The time-dependent forcing protocol `pfunc` describes the evolution of the parameter defined
by `pidx` over a time interval `ptspan`. The parameter `p` that is swept is defined by
`current_parameters(ds)[pidx]`.

The forcing protocol can be specified directly via a `[RateFunction](@def)` referred to as
`rf = RateFunction(pfunc, ptspan)`. The forcing protocol can be adjusted by using the
keyword arguments `forcing_start_time`, `forcing_length`, and `forcing_scale` when
constructing the `RateSystem`. This allows to repurpose the same forcing protocol `rf`.

## Keyword arguments
- `forcing_start_time = ptspan[1]`: Time when parameter shift starts
  (before this, the resulting system will be autonomous)
- `forcing_length = ptspan[2] - ptspan[1]`: Time-interval over which `pfunc(ptspan)` is
  spread out (for window_length > ptspan[2] - ptspan[1]) or squeezed into
  (for window_length < ptspan[2] - ptspan[1])
- `forcing_scale = 1.0`: Amplitude of the protocol.
- `t0 = 0.0`: Initial time of the resulting non-autonomous system

"""

struct RateSystem{S,P,F,T} <: ContinuousTimeDynamicalSystem
    "Underlying continuous-time dynamical system"
    system::S
    "The index of the parameter being forced"
    pidx::P
    "Time-dependent parametric forcing protocol"
    rate_function::RateFunction{F,T}
    "Time when parameter shift starts"
    forcing_start_time::T
    "Time-interval of the protocol"
    forcing_length::T
    "Amplitude of the protocol"
    forcing_scale::T

    function RateSystem(
        ds::ContinuousTimeDynamicalSystem,
        rf::RateFunction{F,T},
        pidx::P;
        forcing_start_time=rf.ptspan[1],
        forcing_length=(rf.ptspan[2] - rf.ptspan[1]),
        forcing_scale=1.0,
        t0=0.0,
    ) where {P,F,T}
        @assert (forcing_start_time >= t0) "The forcing cannot start before the system start time t0, but your forcing_start_time ($(forcing_start_time)) < t0 ($(t0))."

        nonauto_ds = apply_ramping(
            ds, rf, pidx, forcing_start_time, forcing_length, forcing_scale, t0
        )

        return new{typeof(nonauto_ds),P,F,T}(
            nonauto_ds, pidx, rf, forcing_start_time, forcing_length, forcing_scale
        )
    end

    function RateSystem(
        ds::ContinuousTimeDynamicalSystem,
        pfunc::F,
        ptspan::Tuple,
        pidx;
        forcing_start_time=ptspan[1],
        forcing_length=(ptspan[2] - ptspan[1]),
        forcing_scale=1.0,
        t0=0.0,
    ) where {F}
        forcing_start_time, forcing_length, forcing_scale, _ = promote(
            forcing_start_time, forcing_length, forcing_scale, ptspan[1]
        )
        T = typeof(forcing_scale)
        ptspan = reinterpret(NTuple{2,T}, ptspan)

        rf = RateFunction{F,T}(pfunc, ptspan)

        rs = RateSystem(ds, rf, pidx; forcing_start_time, forcing_length, forcing_scale, t0)
        return rs
    end
end

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
DynamicalSystemsBase.isdeterministic(rs::RateSystem) = false

(rs::RateSystem)(t::Real) = rs.system.integ(t)

"""
Creates a piecewise constant function in alignment with the entries of the RateFunction and
the parameter value of the underlying autonomous system
"""
function p_modified(
    t::Real, rf::RateFunction, p0, forcing_start_time, forcing_length, forcing_scale
)
    p = rf.pfunc
    section_start = rf.ptspan[1]
    section_end = rf.ptspan[2]

    # making the function piecewise constant with range [p0,p0+forcing_scale]
    q = if t ≤ forcing_start_time
        return p0
    else
        if t < forcing_start_time + forcing_length
            # performing the time shift corresponding to stretching/squeezing
            time_shift =
                ((section_end - section_start) / forcing_length) *
                (t - forcing_start_time) + section_start
            return p0 +
                   forcing_scale * (p(time_shift) - p(section_start)) /
                   (p(section_end) - p(section_start))
        else
            return p0 + forcing_scale
        end
    end

    return q
end

"returns a continuous time dynamical system with modified drift field"
function apply_ramping(ds, rc, pidx, forcing_start_time, forcing_length, forcing_scale, t0)
    f = dynamic_rule(ds)
    params = current_parameters(ds)
    p0 = params[pidx]

    # this should be rewritten with the Callback mechanism from SciML
    function f_new(u, p, t)
        pvalue = p_modified(t, rc, p0, forcing_start_time, forcing_length, forcing_scale)
        if p isa Union{AbstractArray,AbstractDict}
            setindex!(p, pvalue, pidx)
        else
            setproperty!(p, pidx, pvalue)
        end
        return f(u, p, t)
    end
    old_prob = deepcopy(referrenced_sciml_prob(ds))
    prob = remake(old_prob; f=f_new, p=params, tspan=(t0, Inf))
    return CoupledODEs(prob, ds.diffeq)
end

function set_forcing_scale!(rs::RateSystem, scale) end

function set_forcing_length!(rs::RateSystem, length) end

function set_forcing_start!(rs::RateSystem, start) end
