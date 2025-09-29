# using Pkg
# Pkg.activate(homedir()*"/PhD/Software/CriticalTransitions.jl/");
# Pkg.develop;

using CriticalTransitions

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


# Now, we want to apply a parameter shift to auto_sys s.t. speed and amplitude of the parameter shift are easily modifiable 
# but the 'shape' of it is always the same - just linearly stretched or squeezed. This can be achieved as follows:

# 1) First specify a section of a function p(t) that you would like to use to ramp a parameter of auto_sys:
p(t) = tanh(t);         # A monotonic function that describes the parameter shift
section_start = -100;   # start of the section of p(t) we want to consider
section_end = 100;      # end   of the section of p(t) we want to consider
rc = RateConfig(p, section_start, section_end)

# 2) Then specify how you would like to use the RateConfig to ramp the pidx'th parameter of auto_sys:
pidx = 1
forcing_start = -50.    # time when the parameter shift should start (before this, the final system will be autonomous)
forcing_length = 105.   # time over which p([section_start, section_end]) is spread out (for window_length > section_end - section_start) or squeezed into (for window_length < section_end - section_start)
forcing_scale = 3.0    # amplitude of the ramping. `p` is then automatically rescaled (it will go from p0 to p0+dp with p0 being the value of the pidx'th parameter of the auto_sys)
t0 = -70.0              # initial time of the resulting non-autonomous system (relevant to later compute trajectories)

# Note: for the given section p([-100, 100]) the critical forcing_length is approximately 100 for p0 = 0 and dp = 3
# Hence, for window_length < 100 the trajectory tips, and for window_length > 100 the trajectory tracks to -4.

nonauto_sys = RateSystem(auto_sys, rc, pidx;
    forcing_start=forcing_start,
    forcing_length=forcing_length,
    forcing_scale=forcing_scale,
    t0=t0)

# Compute trajectories
T = forcing_length + 40.0;        # simulation time
auto_traj = trajectory(auto_sys, T, x0);   
nonauto_traj = trajectory(nonauto_sys.system, T, x0);

# Plot trajectories
using CairoMakie
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
p_auto = [0.];  
auto_sys = CoupledODEs(f,x0,p_auto);

p(t) = tanh(t);         # A monotonic function that describes the parameter shift
section_start = -100;   # start of the section of p(t) we want to consider
section_end = 100;      # end   of the section of p(t) we want to consider
rc = RateConfig(p, section_start, section_end)

# 2) Then specify how you would like to use the RateConfig to ramp the pidx'th parameter of auto_sys:
pidx = 1
forcing_start = -100.    # time when the parameter shift should start (before this, the final system will be autonomous)
forcing_length = 200.   # time over which p([section_start, section_end]) is spread out (for window_length > section_end - section_start) or squeezed into (for window_length < section_end - section_start)
forcing_scale = 1.0    # amplitude of the ramping. `p` is then automatically rescaled (it will go from p0 to p0+dp with p0 being the value of the pidx'th parameter of the auto_sys)
t0 = forcing_start             # initial time of the resulting non-autonomous system (relevant to later compute trajectories)

nonauto_sys = RateSystem(auto_sys, rc, pidx; forcing_start=forcing_start, forcing_length=forcing_length, forcing_scale=forcing_scale, t0=t0)

# Compute trajectory
T = forcing_length + 40.0;        # simulation time
@time nonauto_traj = trajectory(nonauto_sys.system, T, x0);



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




### plot for comparison
fig = Figure(); 
axs = Axis(fig[1,1],xlabel="t",ylabel="x");
lines!(axs,nonauto_traj[2],nonauto_traj[1][:,1],linewidth=2,label=L"\text{RateSyst-Nonautonomous system}");
lines!(axs,nonautoExpl_traj[2],nonautoExpl_traj[1][:,1],linewidth=2,label=L"\text{Expl-Nonautonomous system}");
axislegend(axs,position=:rc,labelsize=10);
fig

