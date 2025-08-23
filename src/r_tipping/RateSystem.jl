# we consider the ODE dxₜ/dt = f(xₜ,p(rt))
# here p = p(t) ∈ Rᵐ is a function containing all the system parameters 

# We ask the user to define: 
#  1) a ContinuousTimeDynamicalSystem that should be investigated and
#  2) a protocol for the time-dependent forcing with the struct RateConfig

# Then we give back the ContinuousTimeDynamicalSystem with the parameter 
# changing according to the rate protocol
"""
    RateConfig

Time-dependent forcing protocol specified by the following fields:
- `p::Function`: forcing function of the form `p(p_parameters, t_start)``
- `p_parameters::Vector`: parameters of the forcing function p(t)
- `r::Float64`: rate parameter
- `t_start::Float64`: start time of protocol
- `t_end::Float64`: end time of protocol

Default values 
==============

t_start = -Inf
t_end = Inf
p_parameters = []
"""
mutable struct RateConfig
    p::Function  # has to be a function that ramps from p(t=-100)=0 to p(t=100)=1
    p_parameters::Vector 
    t_start::Float64
    ramp_t_length::Float64
    dp::Float64
end


# convenience functions

#function RateConfig(p::Function, p_parameters::Vector, r::Float64)
#    RateConfig(p, p_parameters, r, -Inf, Inf)
#end
#RateConfig(p::Function, r::Float64)=RateConfig(p, [], r, -Inf, Inf)

function modified_drift(
    u,
    p_parameters,
    t,
    ds::ContinuousTimeDynamicalSystem,
    p::Function,
    t_start::Float64,
    ramp_t_length::Float64,
    dp::Float64;
    kwargs...,
)

    p̃ = if t ≤ t_start
        dp*p(p_parameters, -100.0; kwargs...)
    elseif t_start < t < t_start+ramp_t_length
        dp*p(p_parameters, -100.0+(200/ramp_t_length)*(t-t_start); kwargs...)
    else
        dp*p(p_parameters, 100.0; kwargs...) 
    end; 
    return dynamic_rule(ds)(u, p̃, t)
end;

"""
    apply_ramping(sys::CoupledODEs, rp::RateConfig, t0=0.0; kwargs...)

Applies a time-dependent [`RateConfig`](@def) to a given autonomous deterministic dynamical system
`sys`, turning it into a non-autonomous dynamical system. The returned [`CoupledODEs`](@ref)
has the explicit parameter time-dependence incorporated.

The returned [`CoupledODEs`](@ref) is autonomous from `t_0` to `t_start`, 
then non-autnonmous from `t_start` to `t_end` with the parameter shift given by the [`RateConfig`](@def),
and again autonomous from `t_end` to the end of the simulation:

`t_0`  autonomous    `t_start`  non-autonomous   `t_end`  autonomous   `∞`

Computing trajectories of the returned [`CoupledODEs`](@ref) can then be done in the same way as for any other [`CoupledODEs`](@ref).
"""
function apply_ramping(auto_sys::CoupledODEs, rp::RateConfig, t0=0.0; kwargs...)
    # we wish to return a continuous time dynamical system with modified drift field

    f(u, p_parameters, t) = modified_drift(
        u, p_parameters, t, auto_sys, rp.p, rp.t_start, rp.ramp_t_length, rp.dp; kwargs...
    )
    prob = remake(referrenced_sciml_prob(auto_sys); f, p=rp.p_parameters, tspan=(t0, Inf))
    nonauto_sys = CoupledODEs(prob, auto_sys.diffeq)
    return nonauto_sys
end


