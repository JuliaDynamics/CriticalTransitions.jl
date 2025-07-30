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
with `λ(t) ∈ Rᵐ` constant before time `t_start` and also after time `t_end`. 

RateProtocol contains the following fields:
- `λ::Function`: forcing function of the form `λ(p, t_start; kwargs...)``
- `p_lambda::Vector`: parameters of forcing function
- `r::Float64`: rate parameter
- `t_start::Float64`: start time of protocol
- `t_end::Float64`: end time of protocol

Default values 
==============

If `t_start` and `t_end` are not provided, they are set to `t_start=-Inf`, and `t_end=Inf`.
If `p_lambda` is not provided, it is set to an empty array `[]`.
"""
mutable struct RateProtocol
    λ::Function
    p_lambda::Vector
    r::Float64
    t_start::Float64
    t_end::Float64
end

# convenience functions

function RateProtocol(λ::Function, p_lambda::Vector, r::Float64)
    RateProtocol(λ, p_lambda, r, -Inf, Inf)
end
RateProtocol(λ::Function, r::Float64)=RateProtocol(λ, [], r, -Inf, Inf)
#RateProtocol(λ::Function,p_lambda::Vector,r::Float64,t_start::Float64)=RateProtocol(λ,p_lambda,r,t_start,Inf)
#RateProtocol(λ::Function,r::Float64,t_start::Float64)=RateProtocol(λ,[],r,t_start,Inf)	

function modified_drift(
    u,
    p,
    t,
    ds::ContinuousTimeDynamicalSystem,
    λ::Function,
    t_start::Float64,
    t_end::Float64,
    r::Float64;
    kwargs...,
)
    if t_start > t_end
        error("Please ensure that t_start ≤ t_end.")
    end

    p̃ = if r*t ≤ t_start
        λ(p, t_start; kwargs...)
    elseif t_start < r*t < t_end
        λ(p, r*t; kwargs...) # the value(s) of λ(rt)
    else
        λ(p, t_end; kwargs...) # the value(s) of λ(rt)
    end; # the value(s) of λ(rt)
    return ds.integ.f(u, p̃, t)
end;

"""
    apply_ramping(sys::CoupledODEs, rp::RateProtocol, t0=0.0; kwargs...)

Applies a time-dependent [`RateProtocol`](@def) to a given autonomous deterministic dynamical system `sys`, 
returning a non-autonomous dynamical system of type [`CoupledODEs`](@ref).

The [`RateProtocol`](@def) replaces the parameters of `sys` by the function `λ(rt)` within the 
time interval `[t_start, t_end]`. Thus, the returned [`CoupledODEs`](@ref) has the explicit parameter time-dependence incorporated and is 
autonomous from `t_0` to `t_start`, non-autnonmous from `t_start` to `t_end` with the parameter shift given by the [`RateProtocol`](@def),
and autonomous from `t_end` to the end of the simulation:

`t_0`  autonomous    `t_start`  non-autonomous   `t_end`  autonomous   `∞`
"""
function apply_ramping(auto_sys::CoupledODEs, rp::RateProtocol, t0=0.0; kwargs...)
    # we wish to return a continuous time dynamical system with modified drift field

    f(u, p, t) = modified_drift(
        u, p, t, auto_sys, rp.λ, rp.t_start, rp.t_end, rp.r; kwargs...
    )
    prob = remake(auto_sys.integ.sol.prob; f, p=rp.p_lambda, tspan=(t0, Inf))
    nonauto_sys = CoupledODEs(prob, auto_sys.diffeq)
    return nonauto_sys
end
