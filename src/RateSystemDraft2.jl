using GLMakie
using CriticalTransitions
using DifferentialEquations

# we consider the ODE dxₜ/dt = f(xₜ,λ(rt))
# here λ = λ(t) ∈ Rᵐ is a function containing all the system parameters 

# We ask the user to define: 
#  1) a system that should be investigated and
#  2) a protocol for λ(rt) by calling the function rate_protocol 


##############################################################################################################
# we write behind the scenes
##############################################################################################################

mutable struct RateProtocol 
	λ::Function
	p_lambda::Vector
	r::Float64
	t_start::Float64
    t_end::Float64
end

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

##############################################################################################################
# user writes (pseudo example - see the bottom of the page for an explicit example)
##############################################################################################################


##############################################################################################################
# test (prototypical model used for studying R-tipping, with critical rate r = 4/3)
##############################################################################################################


# autonomous system

function f(u,p,t)
    x = u[1]
    λ = p[1]
    dx = (x+λ)^2 - 1
    return SVector{1}(dx)
end;

lambda = 0.; 
p = [lambda];
 
z₀ = [-1.];

Δt = 1.e-3;
diffeq = (alg = AutoVern9(Rodas5(autodiff=false)),abstol=1.e-16,reltol=1.e-16);

auto_sys = CoupledODEs(f,z₀,p;diffeq)

#  time-dependent parameter function

function λ(p,t)
    λ_max = p[1]
    lambda = (λ_max/2)*(tanh(λ_max*t/2)+1)
    return SVector{1}(lambda)
end;

# calling RateSystem

r = 4/3-0.02; 
λ_max = 3.; p_lambda = [λ_max]; 

#t_start = -Inf; t_end = Inf;
t0 = -10.;

## test v2

rp = RateProtocol(λ,p_lambda,r);

nonauto_sys = RateSystem(auto_sys,rp,t0);

#nonauto_sys = RateSystem(auto_sys,λ,r,t0;p_lambda,t_start,t_end);

T = 20.;

auto_traj = trajectory(auto_sys,T,z₀);
nonauto_traj = trajectory(nonauto_sys,T,z₀);

fig = Figure(); axs = Axis(fig[1,1]);
lines!(axs,t0.+auto_traj[2],auto_traj[1][:,1],linewidth=2,label=L"\text{Autonomous system}");
lines!(axs,nonauto_traj[2],nonauto_traj[1][:,1],linewidth=2,label=L"\text{Nonautonomous system}");
axislegend(axs,position=:rc,labelsize=40);
ylims!(-4.5,1.5);
fig