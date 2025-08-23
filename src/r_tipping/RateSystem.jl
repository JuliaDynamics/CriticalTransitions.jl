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
- `p::Function`: forcing function of the form `p(p, t_start; kwargs...)``
- `p_parameters::Vector`: parameters of the forcing function p(t)
- `r::Float64`: rate parameter
- `t_start::Float64`: start time of protocol
- `t_end::Float64`: end time of protocol

Default values 
==============

t_start=-Inf
t_end=Inf
p_parameters = []
"""
mutable struct RateConfig
    p::Function
    p_parameters::Vector
    r::Float64
    t_start::Float64
    t_end::Float64
end

# convenience functions

function RateConfig(p::Function, p_parameters::Vector, r::Float64)
    RateConfig(p, p_parameters, r, -Inf, Inf)
end
RateConfig(p::Function, r::Float64)=RateConfig(p, [], r, -Inf, Inf)

function modified_drift(
    u,
    par,
    t,
    ds::ContinuousTimeDynamicalSystem,
    p::Function,
    t_start::Float64,
    t_end::Float64,
    r::Float64;
    kwargs...,
)
    if t_start > t_end
        error("Please ensure that t_start ≤ t_end.")
    end

    p̃ = if r*t ≤ t_start
        p(par, t_start; kwargs...)
    elseif t_start < r*t < t_end
        p(par, r*t; kwargs...) # the value(s) of p(rt)
    else
        p(par, t_end; kwargs...) # the value(s) of p(rt)
    end; # the value(s) of p(rt)
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

    f(u, par, t) = modified_drift(
        u, par, t, auto_sys, rp.p, rp.t_start, rp.t_end, rp.r; kwargs...
    )
    prob = remake(referrenced_sciml_prob(auto_sys); f, p=rp.p_parameters, tspan=(t0, Inf))
    nonauto_sys = CoupledODEs(prob, auto_sys.diffeq)
    return nonauto_sys
end

