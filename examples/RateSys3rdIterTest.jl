using CriticalTransitions
using CairoMakie

function f(u,p,t) # out-of-place
    x = u[1]
    dx = (x+p[1])^2 - 1
    return SVector{1}(dx)
end

x0 = [0.]
auto_sys = CoupledODEs(f,x0,[0.0])

p(t) = tanh(t)

p_plotvals = [p(t)[1] for t in -10.0:0.1:10.0]
figp = Figure(); axsp = Axis(figp[1,1],xlabel="t",ylabel=L"p")
lines!(axsp,-10.0:0.1:10.0,p_plotvals,linewidth=2,label=L"p(t)")
axislegend(axsp,position=:rc,labelsize=10)
figp


pidx=1
section_start = -3    # start time of non-autonomous part
section_end = 3       # start time of non-autonomous part
window_start = 0
window_length = 50
dp=1                    # strength of the paramter ramping

rc = CriticalTransitions.RateConfig(pidx,p,section_start, section_end,window_start, window_length, dp)


t0 = -10.0      # initial time of the system
nonauto_sys = apply_ramping(auto_sys, rc, t0);


T = 100.0        # final simulation time
auto_traj = trajectory(auto_sys, T, x0)
nonauto_traj = trajectory(nonauto_sys, T, x0);


fig = Figure(); axs = Axis(fig[1,1],xlabel="t",ylabel="x")
lines!(axs,t0.+auto_traj[2],auto_traj[1][:,1],linewidth=2,label=L"\text{Autonomous system}")
lines!(axs,nonauto_traj[2],nonauto_traj[1][:,1],linewidth=2,label=L"\text{Nonautonomous system}")
axislegend(axs,position=:rc,labelsize=10)
fig

