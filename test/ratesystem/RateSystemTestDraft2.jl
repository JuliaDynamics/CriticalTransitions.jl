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

# fig = Figure(); axs = Axis(fig[1,1]);
# lines!(axs,t0.+auto_traj[2],auto_traj[1][:,1],linewidth=2,label=L"\text{Autonomous system}");
# lines!(axs,nonauto_traj[2],nonauto_traj[1][:,1],linewidth=2,label=L"\text{Nonautonomous system}");
# axislegend(axs,position=:rc,labelsize=40);
# ylims!(-4.5,1.5);
# fig