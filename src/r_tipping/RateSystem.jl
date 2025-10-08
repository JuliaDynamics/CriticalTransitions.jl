"""
    RateSystem
 
"""
struct RateSystem{IIP,D,I,P,F,T} <: ContinuousTimeDynamicalSystem
    integ::I
    # things we can't recover from `integ`
    p0::P
    diffeq # isn't parameterized because it is only used for display
    forcing::F
end

"""
    RateSystem(sys::ContinuousTimeDynamicalSystem, rc::RateConfig, pidx; kwargs...)

Creates a `RateSystem` type from an autonomous dynamical system `sys` and a time-dependent
parametric forcing protocol of `RateConfig` type.

## Keyword arguments
- `forcing_start = RateConfig.interval[1]`: Time when parameter shift starts (before this, the resulting system will be autonomous)
- `forcing_length = RateConfig.interval[2] - RateConfig.interval[1]`: Time-interval over which RateConfig.pfunc([RateConfig.interval[1], RateConfig.interval[2]]) is spread out (for window_length > RateConfig.interval[2] - RateConfig.interval[1]) or squeezed into (for window_length < RateConfig.interval[2] - RateConfig.interval[1])
- `forcing_scale = 1.0`: Amplitude of the ramping. The ramping is then automatically rescaled 
- `t0 = 0.0`: Initial time of the resulting non-autonomous system (relevant to later compute trajectories)

## Description
The returned `RateSystem` is 
autonomous before `forcing_start`, 
non-autnonmous from `forcing_start` to `forcing_start+forcing_length` with the parameter shift given by the [`RateConfig`](@def), and again 
autonomous after `forcing_start+forcing_length`.
"""
function RateSystem(
    auto_sys::ContinuousTimeDynamicalSystem,
    rc::RateConfig,
    pidx;
    forcing_start=rc.interval[1],
    forcing_length=(rc.interval[2] - rc.interval[1]),
    forcing_scale=1.0,
    t0=0.0,
)
    params = deepcopy(current_parameters(auto_sys))
    p0 = params[pidx]

    @assert (forcing_start >= t0) "The forcing cannot start before the system start time t0, but your forcing_start ($(forcing_start)) < t0 ($(t0))."

    system = apply_ramping(
        auto_sys, rc, pidx, p0, params, forcing_start, forcing_length, forcing_scale, t0
    )

    prob = ODEProblem(system)

    return RateSystem(
        prob, system.diffeq; rc, pidx, forcing_start, forcing_length, forcing_scale, t0
    )
end

function RateSystem(
    prob::ODEProblem,
    diffeq=DEFAULT_DIFFEQ;
    rc=RateConfig(:linear),
    pidx=1,
    forcing_start=rc.interval[1],
    forcing_length=(rc.interval[2] - rc.interval[1]),
    forcing_scale=0.0,
    t0=0.0,
    special_kwargs...,
)
    if haskey(special_kwargs, :diffeq)
        throw(
            ArgumentError(
                "`diffeq` is given as positional argument when an ODEProblem is provided."
            ),
        )
    end
    IIP = isinplace(prob)
    D = length(prob.u0)
    P = typeof(prob.p)
    if prob.tspan === (nothing, nothing)
        # If the problem was made via MTK, it is possible to not have a default timespan.
        U = eltype(prob.u0)
        prob = SciMLBase.remake(prob; tspan=(U(t0), U(Inf)))
    end
    solver, remaining = _decompose_into_solver_and_remaining(diffeq)
    integ = __init(
        prob,
        solver;
        remaining...,
        special_kwargs...,
        # Integrators are used exclusively iteratively. There is no reason to save anything.
        save_start=false,
        save_end=false,
        save_everystep=false,
        # DynamicalSystems.jl operates on integrators and `step!` exclusively,
        # so there is no reason to limit the maximum time evolution
        maxiters=Inf,
    )
    forcing = (pidx, rc, forcing_start, forcing_length, forcing_scale)
    return RateSystem{IIP,D,typeof(integ),P,typeof(forcing)}(
        integ, deepcopy(prob.p), diffeq, forcing
    )
end

function forcing(sys::RateSystem, t)
    return p_modified(t, sys.forcing, prob.p, forcing_start, forcing_length, forcing_scale)
end

## Creates a piecewise constant function in alignment with the entries of the RateConfig and the parameter value of the underlying autonomous system
function p_modified(
    t::Real, rc::RateConfig, p0, forcing_start, forcing_length, forcing_scale
)
    p = rc.pfunc
    section_start = rc.interval[1]
    section_end = rc.interval[2]

    # making the function piecewise constant with range [p0,p0+forcing_scale]    
    q = if t ≤ forcing_start
        return p0
    else
        if t < forcing_start+forcing_length
            # performing the time shift corresponding to stretching/squeezing 
            time_shift =
                ((section_end-section_start)/forcing_length)*(t-forcing_start)+section_start
            return p0 +
                   forcing_scale*(
                p(time_shift)-p(section_start)
            )/(p(section_end)-p(section_start))
        else
            return p0 + forcing_scale
        end
    end

    return q
end

# returns a continuous time dynamical system with modified drift field
function apply_ramping(
    auto_sys, rc, pidx, p0, params, forcing_start, forcing_length, forcing_scale, t0
)
    f = deepcopy(dynamic_rule(auto_sys))

    function f_new(u, p, t)
        pvalue = p_modified(t, rc, p0, forcing_start, forcing_length, forcing_scale)
        if p isa Union{AbstractArray,AbstractDict}
            setindex!(p, pvalue, pidx)
        else
            setproperty!(p, pidx, pvalue)
        end
        return f(u, p, t)
    end

    prob = remake(
        deepcopy(referrenced_sciml_prob(auto_sys)); f=f_new, p=params, tspan=(t0, Inf)
    )
    nonauto_sys = CoupledODEs(prob, auto_sys.diffeq)
    return nonauto_sys
end

function set_forcing_scale!(rs::RateSystem, scale) end

function set_forcing_length!(rs::RateSystem, length) end

function set_forcing_start!(rs::RateSystem, start) end
