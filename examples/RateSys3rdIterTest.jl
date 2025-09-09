# using Pkg
# Pkg.activate(homedir()*"/PhD/Software/CriticalTransitions.jl/");
# Pkg.develop;

using CriticalTransitions
using CairoMakie


# First define an autonomous system of which we later want to make one parameter time-dependent
function f(u,p,t) # out-of-place
    x = u[1]
    λ = p[1]
    dx = (x+λ)^2 - 1
    return SVector{1}(dx)
end;

x0 = [-1.0];
p_auto = [0.];  
auto_sys = CoupledODEs(f,x0,p_auto);


# Set up the parameter ramping:
# 1) First specify 
#       - the index pidx of the parameter of the autonomous system that you would like to apply the shift to
# 2) Then the shape of the parameter shift you would like to consider by giving
#       - a monotonic function p(t) describing the shape of the shift
#       - an interval [section_start, section_end] over which p(t) should be considered
#
# 3) Specify how you would like to use this section p([section_start, section_end]) to shift the pidx'th parameter 
#    of the autonomous system. This is done by specifying:
#       - a window_start, i.e. the time when the parameter shift should start (before this, the system is autonomous)
#       - a window_length over which p([section_start, section_end]) is 
#           spread out (for window_length > section_end - section_start) or 
#           squeezed into (for window_length < section_end - section_start)
#       - an amplitude dp. Then p(t) is automatically scaled s.t. p(section_end) - p(section_start) = dp

pidx=1;                 # Index of the parameter that is shifted 
p(t) = tanh(t);         # Function that describes the parameter shift

# plot p(t) to see which part of it we would like to consider
p_plotvals = [p(t) for t in -10.0:0.1:10.0];
figp = Figure(); 
axsp = Axis(figp[1,1],xlabel="t",ylabel=L"p");
lines!(axsp,-10.0:0.1:10.0,p_plotvals,linewidth=2,label=L"p(t)");
axislegend(axsp,position=:rc,labelsize=10);
figp

# now specify the section of p(t) you would like to consider
section_start = -100;    # start time of the section of the tanh we want to consider
section_end = 100;       # start time of the section of the tanh we want to consider

window_start = -50.;    # time when the parameter shift should start (before this, the system is autonomous)

# For the section p([-100, 100]) the critical window_length is approximately 100 for p0 = 0 and dp = 3
# Hence for window_length < 100, the trajectory tips
# and   for window_length > 100, the trajectory tracks to -4
# Thus, choose one of the following 2 lines:
# window_length = 95; # tips!
window_length = 105; # tracks!

dp=3;         # amplitude of the parameter ramping (i.e. it will go from p0 to p0+dp)

# set up the nonauto_sys by combining all this information in a RateConfig
rc = CriticalTransitions.RateConfig(pidx,p,section_start, section_end, window_start, window_length, dp);
t0 = window_start - 20.0;      # initial time of the system
nonauto_sys = apply_ramping(auto_sys, rc, t0);

# Compute trajectories
T = window_length + 40.0;        # simulation time
auto_traj = trajectory(auto_sys, T, x0);   
nonauto_traj = trajectory(nonauto_sys, T, x0);

# Plot trajectories
fig = Figure(); 
axs = Axis(fig[1,1],xlabel="t",ylabel="x");
lines!(axs,t0.+auto_traj[2],auto_traj[1][:,1],linewidth=2,label=L"\text{Autonomous system}");
lines!(axs,nonauto_traj[2],nonauto_traj[1][:,1],linewidth=2,label=L"\text{Nonautonomous system}");
axislegend(axs,position=:rc,labelsize=10);
fig










####################
### BENCHMARKING ###
####################

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

pidx=1;
p(t) = tanh(t);
section_start = -100;   # start time of non-autonomous part
section_end = 100;      # start time of non-autonomous part
window_start = -100.;
window_length = 200;    # tracks!
dp=1;                   # strength of the paramter ramping
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



### plot for comparison
fig = Figure(); 
axs = Axis(fig[1,1],xlabel="t",ylabel="x");
lines!(axs,nonauto_traj[2],nonauto_traj[1][:,1],linewidth=2,label=L"\text{RateSyst-Nonautonomous system}");
lines!(axs,nonautoExpl_traj[2],nonautoExpl_traj[1][:,1],linewidth=2,label=L"\text{Expl-Nonautonomous system}");
axislegend(axs,position=:rc,labelsize=10);
fig

