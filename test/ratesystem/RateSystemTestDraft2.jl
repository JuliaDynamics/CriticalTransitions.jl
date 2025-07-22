##############################################################################################################
# test (prototypical model used for studying R-tipping, with critical rate r = 4/3)
##############################################################################################################

# autonomous system
function f(u,p,t)
    x = u[1]
    λ = p[1]
    dx = (x+λ)^2 - 1
    return SVector{1}(dx)
end

lambda = 0.
p = [lambda]
x0 = [-1.]

Δt = 1.e-3
#diffeq = (alg = AutoVern9(Rodas5(autodiff=false)),abstol=1.e-16,reltol=1.e-16)
#auto_sys = CoupledODEs(f,x0,p;diffeq)
auto_sys = CoupledODEs(f,x0,p)


# now we make the system non-autonomous

#  We define a time-dependent parameter function
function λ(p,t)
    λ_max = p[1]
    lambda = (λ_max/2)*(tanh(λ_max*t/2)+1)
    return SVector{1}(lambda)
end

# we set the parameters for RateProtocol
λ_max = 3.
p_lambda = [λ_max]
r = 4/3-0.02 # r just below critical rate
t_start = -Inf # start time of non-autonomous part
t_end = Inf    # end time of non-autonomous part
# And the initial time of the system
t0 = -10.

rp = CriticalTransitions.RateProtocol(λ,p_lambda,r,t_start,t_end)
nonauto_sys = RateSystem(auto_sys,rp,t0)

T = 20. # final simulation time
auto_traj = trajectory(auto_sys,T,x0)
nonauto_traj = trajectory(nonauto_sys,T,x0)

fig = Figure(); axs = Axis(fig[1,1])
lines!(axs,t0.+auto_traj[2],auto_traj[1][:,1],linewidth=2,label=L"\text{Autonomous system}")
lines!(axs,nonauto_traj[2],nonauto_traj[1][:,1],linewidth=2,label=L"\text{Nonautonomous system}")
axislegend(axs,position=:rc,labelsize=10)
fig