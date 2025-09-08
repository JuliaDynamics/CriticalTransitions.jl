# using Pkg
# Pkg.activate(homedir()*"/PhD/Software/CriticalTransitions.jl/");
# Pkg.develop;

using CriticalTransitions
using CairoMakie

function f(u,p,t) # out-of-place
    x = u[1]
    λ = p[1]
    dx = (x+λ)^2 - 1
    return SVector{1}(dx)
end;

x0 = [-1.0];
p_auto = [0.];  # p_auto somehow needs to be redefined every time.
auto_sys = CoupledODEs(f,x0,p_auto);

p(t) = tanh(t);

p_plotvals = [p(t) for t in -10.0:0.1:10.0];
figp = Figure(); 
axsp = Axis(figp[1,1],xlabel="t",ylabel=L"p");
lines!(axsp,-10.0:0.1:10.0,p_plotvals,linewidth=2,label=L"p(t)");
axislegend(axsp,position=:rc,labelsize=10);
figp

pidx=1;
section_start = -100;    # start time of non-autonomous part
section_end = 100;       # start time of non-autonomous part
dp=3;                    # strength of the paramter ramping

# for this nearly full section of the tanh function the critical window_length is approximately 100 for p0 = 0 and dp = 3 :)
# hence when choosing a window_length (roughly) less than 100, the pullback attractor should tip in forward time to +∞
# and when choosing a window_length (roughly) greater than 100, the pullback attractor should track in forward time to -4

#window_start = -50.;
#window_length = 95; # tips!

window_start = -50.;
window_length = 105; # tracks!

rc = CriticalTransitions.RateConfig(pidx,p,section_start, section_end, window_start, window_length, dp);

t0 = window_start - 20.0;      # initial time of the system
nonauto_sys = apply_ramping(auto_sys, rc, t0);

T = window_length + 40.0;        # total simulation time
#auto_traj = trajectory(auto_sys, T, x0);
nonauto_traj = trajectory(nonauto_sys, T, x0);

fig = Figure(); 
axs = Axis(fig[1,1],xlabel="t",ylabel="x");
#lines!(axs,t0.+auto_traj[2],auto_traj[1][:,1],linewidth=2,label=L"\text{Autonomous system}");
lines!(axs,nonauto_traj[2],nonauto_traj[1][:,1],linewidth=2,label=L"\text{Nonautonomous system}");
axislegend(axs,position=:rc,labelsize=10);
#ylims!(axs,-2.0,10.0)
fig






### BENCHMARKING
## Version 1) current RateSyst setting
function f(u,p,t) # out-of-place
    x = u[1]
    λ = p[1]
    dx = (x+λ)^2 - 1
    return SVector{1}(dx)
end;

x0 = [-1.0];
p_auto = [0.];  # p_auto somehow needs to be redefined every time.
auto_sys = CoupledODEs(f,x0,p_auto);

p(t) = tanh(t);

p_plotvals = [p(t) for t in -10.0:0.1:10.0];
figp = Figure(); 
axsp = Axis(figp[1,1],xlabel="t",ylabel=L"p");
lines!(axsp,-10.0:0.1:10.0,p_plotvals,linewidth=2,label=L"p(t)");
axislegend(axsp,position=:rc,labelsize=10);
figp

pidx=1;
section_start = -100;    # start time of non-autonomous part
section_end = 100;       # start time of non-autonomous part
dp=1;                    # strength of the paramter ramping

window_start = -100.;
window_length = 200; # tracks!

rc = CriticalTransitions.RateConfig(pidx,p,section_start, section_end, window_start, window_length, dp);

t0 = window_start;      # initial time of the system
nonauto_sys = apply_ramping(auto_sys, rc, t0);

T = window_length;        # total simulation time
#auto_traj = trajectory(auto_sys, T, x0);
@time nonauto_traj = trajectory(nonauto_sys, T, x0)



## Version 2) hardcoded

function fexpl(u,p,t) # out-of-place
    x = u[1]
    λ = p[1]

    dx = (x+(tanh(t)+1)/2)^2 - 1
    return SVector{1}(dx)
end;

x0 = [-1.0];
p_auto = [0.];  # p_auto somehow needs to be redefined every time.
t0 = -100
nonauto_sysexpl = CoupledODEs(fexpl,x0,p_auto, t0=t0);

@time nonautoExpl_traj = trajectory(nonauto_sysexpl, T, x0)

fig = Figure(); 
axs = Axis(fig[1,1],xlabel="t",ylabel="x");
#lines!(axs,t0.+auto_traj[2],auto_traj[1][:,1],linewidth=2,label=L"\text{Autonomous system}");
lines!(axs,nonauto_traj[2],nonauto_traj[1][:,1],linewidth=2,label=L"\text{RateSyst-Nonautonomous system}");
lines!(axs,nonautoExpl_traj[2],nonautoExpl_traj[1][:,1],linewidth=2,label=L"\text{Expl-Nonautonomous system}");
axislegend(axs,position=:rc,labelsize=10);
#ylims!(axs,-2.0,10.0)
fig

