using CriticalTransitions

# First define an autonomous system of which we later want to make one parameter time-dependent
function f(u, p, t) # out-of-place
    x = u[1]
    λ = p[1]
    dx = (x+λ)^2 - 1
    return SVector{1}(dx)
end;

x0 = [-1.0];
p_auto = [0.0];
auto_sys = CoupledODEs(f, x0, p_auto);

# Now, we want to apply a parameter shift to auto_sys s.t. speed and amplitude of the parameter shift are easily modifiable
# but the 'shape' of it is always the same - just linearly stretched or squeezed. This can be achieved as follows:

# 1) First specify a section of a function p(t) that you would like to use to ramp a parameter of auto_sys:
p(t) = tanh(t);         # A monotonic function that describes the parameter shift
interval = (-5, 5)
rc = RateConfig(p, interval)

# 2) Then specify how you would like to use the RateConfig to ramp the pidx'th parameter of auto_sys:
pidx = 1                # Index of the parameter-container of auto_sys that you want to ramp
forcing_start = -50.0    # Time when the parameter shift should start (before this, the final system will be autonomous)
forcing_length = 105.0   # Time-interval over which p(interval) is spread out or squeezed
forcing_scale = 3.0     # Amplitude of the ramping. `p` is then automatically rescaled
t0 = -70.0              # Initial time of the resulting non-autonomous system (relevant to later compute trajectories)

RateSys = RateSystem(
    auto_sys,
    rc,
    pidx;
    forcing_start=forcing_start,
    forcing_length=forcing_length,
    forcing_scale=forcing_scale,
    t0=t0,
)

# Compute trajectories
T = forcing_length + 40.0;                      # length of the trajectory that we want to compute
auto_traj = trajectory(auto_sys, T, x0);
nonauto_traj = trajectory(RateSys.system, T, x0);

# Plot trajectories
using CairoMakie
fig = Figure();
axs = Axis(fig[1, 1]; xlabel="t", ylabel="x");
lines!(
    axs,
    t0 .+ auto_traj[2],
    auto_traj[1][:, 1];
    linewidth=2,
    label=L"\text{Autonomous system}",
);
lines!(
    axs,
    nonauto_traj[2],
    nonauto_traj[1][:, 1];
    linewidth=2,
    label=L"\text{Nonautonomous system}",
);
axislegend(axs; position=:rc, labelsize=10);
fig

# Plot ramped parameter
figPar = Figure();
axsPar = Axis(figPar[1, 1]; xlabel="t", ylabel="p(t)");
lines!(
    axsPar,
    range(t0, t0+T; length=Int(T+1)),
    RateSys.forcing.(range(t0, t0+T; length=Int(T+1)));
    linewidth=2,
);
figPar

# Note: for the given section p([-100, 100]) the critical forcing_length is approximately 100 for p0 = 0 and forcing_scale = 3
# Hence, for window_length < 100 the trajectory tips, and for window_length > 100 the trajectory tracks to -4.

# %%
using BenchmarkTools
####################
### BENCHMARKING ###
####################

## Version 1) current RateSyst setting
function f(u, p, t) # out-of-place
    x = u[1]
    λ = p[1]
    dx = (x+λ)^2 - 1
    return SVector{1}(dx)
end;
x0 = [-1.0];
p_auto = [0.0];
auto_sys = CoupledODEs(f, x0, p_auto);

p(t) = tanh(t);         # A monotonic function that describes the parameter shift
interval = (-100,100)
rc = RateConfig(p, interval)

# 2) Then specify how you would like to use the RateConfig to ramp the pidx'th parameter of auto_sys:
pidx = 1
forcing_start = -100.0    # time when the parameter shift should start (before this, the final system will be autonomous)
forcing_length = 200.0   # time over which p(interval) is spread out or squeezed
forcing_scale = 1.0    # amplitude of the ramping. `p` is then automatically rescaled (it will go from p0 to p0+dp with p0 being the value of the pidx'th parameter of the auto_sys)
t0 = forcing_start             # initial time of the resulting non-autonomous system (relevant to later compute trajectories)

nonauto_sys = RateSystem(
    auto_sys,
    rc,
    pidx;
    forcing_start=forcing_start,
    forcing_length=forcing_length,
    forcing_scale=forcing_scale,
    t0=t0,
)

# Compute trajectory
T = forcing_length + 40.0;        # simulation time
@btime nonauto_traj = trajectory($(nonauto_sys.system), T, x0);

## Version 2) hardcoded
function fexpl(u, p, t) # out-of-place
    x = u[1]
    λ = p[1]
    dx = (x+(tanh(t)+1)/2)^2 - 1
    return SVector{1}(dx)
end;
x0 = [-1.0];
p_auto = [0.0];  # p_auto somehow needs to be redefined every time.
t0 = -100
nonauto_sysexpl = CoupledODEs(fexpl, x0, p_auto; t0=t0);
@btime nonautoExpl_traj = trajectory($nonauto_sysexpl, T, x0);

### plot for comparison
fig = Figure();
axs = Axis(fig[1, 1]; xlabel="t", ylabel="x");
lines!(
    axs,
    nonauto_traj[2],
    nonauto_traj[1][:, 1];
    linewidth=2,
    label=L"\text{RateSyst-Nonautonomous system}",
);
lines!(
    axs,
    nonautoExpl_traj[2],
    nonautoExpl_traj[1][:, 1];
    linewidth=2,
    label=L"\text{Expl-Nonautonomous system}",
);
axislegend(axs; position=:rc, labelsize=10);
fig
