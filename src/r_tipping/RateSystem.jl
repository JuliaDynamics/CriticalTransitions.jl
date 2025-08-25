# we consider the ODE dxₜ/dt = f(xₜ,p(rt))
# here p = p(t) ∈ Rᵐ is a function containing all the system parameters 

# We ask the user to define: 
#  1) a ContinuousTimeDynamicalSystem that should be investigated and
#  2) a protocol for the time-dependent forcing with the struct RateConfig

# Then we give back the ContinuousTimeDynamicalSystem with the parameter 
# changing according to the rate protocol
"""
    RateConfig

Time-dependent forcing protocol containing the information to apply a parameter shift to an autonomous system.

Fields
==============

- p: monotonic function which describes the time-dependent parametric forcing 
- p_parameters: the vector of parameters which are associated with p
- t_pstart: the parameter values of the past limit system are given by p(p_parameteters,t_pstart)
- t_pend: the parameter values of the future limit system are given by p(p_parameters,t_pend) 
- t_start: the explicit value of time at which the nonautonomous dynamics starts 
- t_ramp_length: the duration of time of nonautonomous dynamics [i.e. within (t_start,t_start+t_ramp_length)] 
- dp: the difference in parameter values attained across the ramping

Default values 
==============

- t_pstart = -100
- t_pend = 100
- p_parameters = []
- t_start = -t_ramp_length/2
- dp = 1
"""
mutable struct RateConfig
    p::Function  
    p_parameters::Vector
    t_pstart::Float64
    t_pend::Float64  
    t_start::Float64
    t_ramp_length::Float64
    dp::Float64
end

## convenience functions

RateConfig(p::Function, t_ramp_length::Float64) = RateConfig(p, [], -100.0, 100.0, -t_ramp_length/2, t_ramp_length, 1.0)
RateConfig(p::Function, t_ramp_length::Float64, dp::Float64) = RateConfig(p, [], -100.0, 100.0, -t_ramp_length/2, t_ramp_length, dp)

RateConfig(p::Function, p_parameters::Vector, t_ramp_length::Float64) = RateConfig(p, p_parameters, -100.0, 100.0, -t_ramp_length/2, t_ramp_length, 1.0)
RateConfig(p::Function, p_parameters::Vector, t_ramp_length::Float64, dp::Float64) = RateConfig(p, p_parameters, -100.0, 100.0, -t_ramp_length/2, t_ramp_length, dp)

## the following function creates a piecewise constant function with respect to t_pstart and t_pend...
## ...which is written such that p̃(t_pend)-p̃(t_pstart) = dp, for the time-dependent entries (otherwise p̃(t_pend)-p̃(t_pstart) = 0 as p̃(t)≡c)

function p_piecewise_scaled(p::Function,p_parameters::Vector,t_pstart::Float64,t_pend::Float64,dp::Float64,t::Float64)
    if t ≤ t_pstart
        return p(p_parameters,t_pstart)
    else
        np = length(p(p_parameters,t_pstart)) # the number of system parameters 
        differences = p(p_parameters,t_pend) .- p(p_parameters,t_pstart)
        nonautonomous_inds = [ii for ii ∈ 1:np if differences[ii] != 0]
        inds = [ii ∈ nonautonomous_inds ? 1/abs(differences[ii]) : 0 for ii ∈ 1:np] # takes value 1/|p(...,t_pend)-p(...,t_start)| when p(...,t_pend)=/=p(...,t_start); zero otherwise
        if t < t_pend
            return p(p_parameters,t_pstart) .+ dp .* inds .* (p(p_parameters,t) .- p(p_parameters,t_pstart))
        else 
            return p(p_parameters,t_pstart) .+ dp .* inds .* (p(p_parameters,t_pend) .- p(p_parameters,t_pstart))
        end
    end
end

function modified_drift(
    u,
    p_parameters,
    t,
    ds::ContinuousTimeDynamicalSystem,
    p::Function,
    t_pstart::Float64,
    t_pend::Float64,
    t_start::Float64,
    t_ramp_length::Float64,
    dp::Float64;
    kwargs...,
)

    q(t) = p_piecewise_scaled(p,p_parameters,t_pstart,t_pend,dp,t)

    time_shift = ((t_pend-t_pstart)/t_ramp_length)*(t-t_start)+t_pstart # such that [t_start,t_start+t_ramp_length] is shifted into [t_pstart,t_pend]

    q̃ = q(time_shift)

    return dynamic_rule(ds)(u, q̃, t)
end;

"""
    apply_ramping(sys::CoupledODEs, rc::RateConfig, t0=0.0; kwargs...)

Applies a time-dependent [`RateConfig`](@def) to a given autonomous deterministic dynamical system
`sys`, turning it into a non-autonomous dynamical system. The returned [`CoupledODEs`](@ref)
has the explicit parameter time-dependence incorporated.

The returned [`CoupledODEs`](@ref) is autonomous from `t_0` to `t_start`, 
then non-autnonmous from `t_start` to `t_start+t_ramp_length` with the parameter shift given by the [`RateConfig`](@def),
and again autonomous from `t_start+t_ramp_length` to the end of the simulation:

`t_0`  autonomous    `t_start`  non-autonomous   `t_start+t_ramp_length`  autonomous   `∞`

Computing trajectories of the returned [`CoupledODEs`](@ref) can then be done in the same way as for any other [`CoupledODEs`](@ref).
"""
function apply_ramping(auto_sys::CoupledODEs, rc::RateConfig, t0=0.0; kwargs...)
    # we wish to return a continuous time dynamical system with modified drift field

    f(u, p_parameters, t) = modified_drift(
        u, p_parameters, t, auto_sys, rc.p, rc.t_pstart, rc.t_pend, rc.t_start, rc.t_ramp_length, rc.dp; kwargs...
    )
    prob = remake(referrenced_sciml_prob(auto_sys); f, p=rc.p_parameters, tspan=(t0, Inf))
    nonauto_sys = CoupledODEs(prob, auto_sys.diffeq)
    return nonauto_sys
end



