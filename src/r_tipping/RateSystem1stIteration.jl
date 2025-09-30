# we consider the ODE dxₜ/dt = f(xₜ,λ(rt))
# here λ = λ(t) ∈ Rᵐ is a function containing all the system parameters 

# We ask the user to define: 
#  1) a ContinuousTimeDynamicalSystem that should be investigated and
#  2) a protocol for the time-dependent forcing with the struct RateProtocol

# Then we give back the ContinuousTimeDynamicalSystem with the parameter 
# changing according to the rate protocol
"""
    RateProtocol

The RateProtocol contains all the fields (information) that allows to take an ODE of the form 
    `dxₜ/dt = f(xₜ, λ)` 
with `λ ∈ Rᵐ` containing all system parameters, and make it an ODE of the form 
    `dxₜ/dt = f(xₜ, λ(rt))`
with `λ(t) ∈ Rᵐ` constant before time `time_interval[1]` and also after time `time_interval[2]`. 

RateProtocol contains the following fields:
- `λ::Function`: forcing function of the form `λ(p, t_start; kwargs...)``
- `p_lambda::Vector`: parameters of the forcing function
- `r::Float64`: rate parameter
- `time_interval::Tuple`: start and end time of protocol

Default values 
==============

time_interval = (-Inf, Inf)
p_lambda = []
"""
mutable struct RateProtocol
    λ::Function
    p_lambda::Vector
    r::Float64
    time_interval::Tuple
end

# convenience functions

function RateProtocol(λ::Function, p_lambda::Vector, r::Float64)
    RateProtocol(λ, p_lambda, r, (-Inf, Inf))
end
RateProtocol(λ::Function, r::Float64)=RateProtocol(λ, [], r, (-Inf, Inf))
#RateProtocol(λ::Function,p_lambda::Vector,r::Float64,t_start::Float64)=RateProtocol(λ,p_lambda,r,t_start,Inf)
#RateProtocol(λ::Function,r::Float64,t_start::Float64)=RateProtocol(λ,[],r,t_start,Inf)	

function modified_drift(
    u,
    p,
    t,
    ds::ContinuousTimeDynamicalSystem,
    λ::Function,
    time_interval::Tuple,
    r::Float64;
    kwargs...,
)
    if time_interval[1] > time_interval[2]
        error("Please ensure that time_interval[1] ≤ time_interval[2].")
    end

    p̃ = if r*t ≤ time_interval[1]
        λ(p, time_interval[1]; kwargs...)
    elseif time_interval[1] < r*t < time_interval[2]
        λ(p, r*t; kwargs...) # the value(s) of λ(rt)
    else
        λ(p, time_interval[2]; kwargs...) # the value(s) of λ(rt)
    end; # the value(s) of λ(rt)
    return dynamic_rule(ds)(u, p̃, t)
end;

"""
    apply_ramping(sys::CoupledODEs, rp::RateProtocol, t0=0.0; kwargs...)

Applies a time-dependent [`RateProtocol`](@def) to a given autonomous deterministic dynamical system `sys`, 
returning a non-autonomous dynamical system of type [`CoupledODEs`](@ref).

The [`RateProtocol`](@def) replaces the parameters of `sys` by the function `λ(rt)` within the 
time interval `time_interval`. Thus, the returned [`CoupledODEs`](@ref) has the explicit parameter time-dependence incorporated and is 
autonomous from `t_0` to `time_interval[1]`, non-autnonmous from `time_interval[1]` to `time_interval[2]` with the parameter shift given by the [`RateProtocol`](@def),
and autonomous from `time_interval[2]` to the end of the simulation:

`t_0`  autonomous    `time_interval[1]`  non-autonomous   `time_interval[2]`  autonomous   `∞`

Computing trajectories of the returned [`CoupledODEs`](@ref) can then be done in the same way as for any other [`CoupledODEs`](@ref).
"""
function apply_ramping(auto_sys::CoupledODEs, rp::RateProtocol, t0=0.0; kwargs...)
    # we wish to return a continuous time dynamical system with modified drift field

    f(u, p, t) = modified_drift(u, p, t, auto_sys, rp.λ, rp.time_interval, rp.r; kwargs...)
    prob = remake(referrenced_sciml_prob(auto_sys); f, p=rp.p_lambda, tspan=(t0, Inf))
    nonauto_sys = CoupledODEs(prob, auto_sys.diffeq)
    return nonauto_sys
end
