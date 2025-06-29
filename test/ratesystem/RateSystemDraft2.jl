# we consider the ODE dxₜ/dt = f(xₜ,λ(rt))
# here λ = λ(t) ∈ Rᵐ is a function containing all the system parameters 

# We ask the user to define: 
#  1) a ContinuousTimeDynamicalSystem that should be investigated and
#  2) a protocol for the time-dependent forcing with the struct RateProtocol

# Then we give back the ContinuousTimeDynamicalSystem with the parameter 
# changing according to the rate protocol
mutable struct RateProtocol 
	λ::Function
	p_lambda::Vector
	r::Float64
	t_start::Float64
    t_end::Float64
end

# convenience functions

RateProtocol(λ::Function,p_lambda::Vector,r::Float64)=RateProtocol(λ,p_lambda,r,-Inf,Inf)
RateProtocol(λ::Function,r::Float64)=RateProtocol(λ,[],r,-Inf,Inf)
#RateProtocol(λ::Function,p_lambda::Vector,r::Float64,t_start::Float64)=RateProtocol(λ,p_lambda,r,t_start,Inf)
#RateProtocol(λ::Function,r::Float64,t_start::Float64)=RateProtocol(λ,[],r,t_start,Inf)	

function modified_drift(u,p,t,ds::ContinuousTimeDynamicalSystem,λ::Function,t_start::Float64,t_end::Float64,r::Float64;
    kwargs...)
	
    if t_start > t_end 
        error("Please ensure that t_start ≤ t_end.")
    end

    p̃ = r*t ≤ t_start ? λ(p,t_start;kwargs...) : t_start < r*t < t_end ? λ(p,r*t;kwargs...) : λ(p,t_end;kwargs...); # the value(s) of λ(rt)
    return ds.integ.f(u,p̃,t)
end;

function RateSystem(auto_sys::ContinuousTimeDynamicalSystem, rp::RateProtocol, t0::Float64;
    kwargs...)
    # we wish to return a continuous time dynamical system with modified drift field

    f(u,p,t) = modified_drift(u,p,t,auto_sys,rp.λ,rp.t_start,rp.t_end,rp.r;kwargs...); 
    prob = remake(auto_sys.integ.sol.prob;f,p=rp.p_lambda,tspan=(t0,Inf));
    nonauto_sys = CoupledODEs(prob,auto_sys.diffeq); 
    return nonauto_sys
end;
