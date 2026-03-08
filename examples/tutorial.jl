# # CriticalTransitions.jl Tutorial

# !!! note "DynamicalSystems.jl and Attractors.jl background recommended"
#     CriticalTransitions.jl is an advanced software for the analysis of critical trasitions
#     in dynamical systems. Due to its advanced nature it is recommended that you have basic
#     familiarity with the DynamicalSystems.jl and Attractors.jl packages, by going through
#     their main tutorials.

# The general workflow of CriticalTransitions.jl consists of two steps, similar to DynamicalSystems.jl:

# 1. Define your specific dynamical system type
#    (either [`RateSystem`](@ref) or [`RandomSystem`](@ref), see below in this tutorial).
# 2. Investigate the system by calling existing functions on it
#    (see [API](@ref), this tutorial, and the Examples entries).

# The picture below showcases the two main routes one can go: rate or random,
# as well as how they integrate with the broader DynamicalSystems.jl library
# and the functionality of CriticalTransitions.jl.

# TODO:.

# ## Types of systems

# There are two main system types supported by this package, both being non-autonomous.
# The first type is systems that are driven by noise, primarily stochastic ordinary differential equations.
# The second type is systems whose parameters change with time deterministically according.

# In both cases one often starts with an autonomous deterministic system,
# which is created following the standard **DynamicalSystems.jl** approach.
# For the scope of this tutorial, this will be the FitzHugh-Nagumo model:

# ```math
# \begin{aligned}
# \frac{du}{dt} &= \frac{1}{\epsilon} \left( -\alpha u^3 + \gamma u - \kappa v + I \right) \\
# \frac{dv}{dt} &= -\beta v + u \ ,
# \end{aligned}
# ```

using CriticalTransitions # re-exports `DynamicalSystemsBase`
using CairoMakie
import Random # hide
Random.seed!(1) # hide

function fitzhugh_nagumo(u,p,t)
    u, v = u
    ϵ, β, α, γ, κ, I = p
    du = (-α*u^3 + γ*u - κ*v + I)/ϵ
    dv = -β*v + u
    return SVector(du, dv)
end

p = [1., 3., 1., 1., 1., 0.] # Parameters (ϵ, β, α, γ, κ, I)
u = [0.1, 0.1]
ds = CoupledODEs(fitzhugh_nagumo, u, p)

# ## RateSystem: creation

# Transforming a deterministic `DynamicalSystem` to a `RateSystem` is straightforward.
# All we have to do is define a forcing profile that dictates how the parameter(s)
# may change with time. The simplest case is a linear ramp, captured by:

fp = ForcingProfile(:linear)

# The forcing profile by itself doesn't say when the forcing starts or ends,
# or how much the parameter will increase overall. It only captures the form
# of the time variability.
# The remaining information is encoded when creating the `RateSystem`:

pidx = 6 # which parameter changes
rs = RateSystem(ds, fp, pidx;
    forcing_start_time = 10,
    forcing_duration = 10,
    forcing_scale = 5
)

# In the above example the current `I` will linearly ramp from 0 to 5 in the time window 10 to 20.
# Being explicit, this is what's going on, given time `t`:

# - `t < forcing_start_time`: the system is autonomous, with parameters given by
#   the underlying autonomous system.
# - `forcing_start_time < t < forcing_start_time + forcing_duration`: the system is
#   non-autonomous with the parameter change given by the `ForcingProfile`,
#   scaled in magnitude by `forcing_scale`.
# - `t > forcing_start_time + forcing_duration`: the system is again autonomous, with
#   parameters fixed at their values attained at the end of the forcing interval.

# Let's simulate both the autonomous and non-autonomous systems to see the difference:

T = 50.0 # Total time
traj_ds, = trajectory(ds, T; Δt = 0.01)
traj_rs, = trajectory(rs, T; Δt = 0.01)
fig = Figure()
tvec = 0:0.01:T
ax = Axis(fig[1, 1]; ylabel="u")
lines!(ax, tvec, traj_ds[:, 1]; label = "autonomous")
lines!(ax, tvec, traj_rs[:, 1]; label = "rate-forced")
axislegend(ax)
fig

# The function [`current_parameters`](@ref) is overloaded for `RateSystem` and can be
# called with a second input `t` to provide the parameters at time `t`. This way we can
# plot the forcing parameter as well:

ps_of_t = current_parameters.(rs, tvec)
I_of_t = getindex.(ps_of_t, pidx)
ax = Axis(fig[2, 1]; xlabel="time", ylabel="I(t)")
lines!(ax, tvec, I_of_t)
fig

# And what if you want multiple parameters to be forced, each with a different profile?
# Achieving this is straightforward. Provide a dictionary mapping the parameter
# index to the corresponding forcing profile! For example:

profiles = Dict(
    6 => ForcingProfile(:linear), # as before
    2 => ForcingProfile(x -> x^2, (0.0, 2.0)), # quadratic!
)

# Note that all profiles start and stop at the same system time, although nothing
# stops you from making a piecewise profile function that has 0s at some of its starting or
# ending portion(s). In any case, we now make a rate system by giving the profiles:

rs = RateSystem(ds, profiles;
    forcing_start_time = 10,
    forcing_duration = 10,
    forcing_scale = 5
)

# TODO: Remains to be done.

# ## RandomSystem: creation

# At the moment, the only case of a `RandomSystem` available is that of
# [`CoupledSDEs`](@ref). The DynamicalSystemsBase.jl has an extensive documentation
# on all the possible ways to create it. This can range from specifying just a noise strength
# (leading to a white noise added to all variables) all the way to completely customized
# noise and dynamics mixing with multiplicative correlated noise.

# Here we will keep things simple and add additive white noise of given strength
# to all variables:

using StochasticDiffEq # required for `CoupledSDEs`
sds = CoupledSDEs(ds, p; noise_strength = 0.2)

# To define more complicated noise processes than simple additive white noise,
# you can specify a custom *noise function* and *covariance matrix* in the `CoupledSDEs`
# definition. For more info, see [`CoupledSDEs`](@ref).
# In any case we now have our stochastic system, so let's visualize a couple of
# trajectories (as this is a stochastic system each one will have different noise
# by default)

fig = Figure()
ax = Axis(fig[1, 1]; xlabel="time", ylabel="u")
for _ in 1:3
    traj_sds, = trajectory(sds, T; Δt = 0.1)
    lines!(ax, 0:0.1:T, traj_sds[:, 1])
end
fig

# ## Usage with the broader DynamicalSystems.jl.

# Both random and rate systems are part of DynamicalSystems.jl and usable directly with it.
# The way to achieve this is to cast the systems back to their deterministic autonomous forms
# (as this is the type of systems the rest of DynamicalSystems.jl covers).
# For example, we may be interested in finding the basins of attraction of the deterministic system

using Attractors

grid = (range(-2, 2; length = 100), range(-2, 2; length = 100),)
mapper = AttractorsViaRecurrences(ds, grid)
boa = basins_of_attraction(mapper, grid)
fig = heatmap_basins_attractors(boa)

# We can visualize a couple of stochastic trajectories on top of the basins of attraction:

ax = content(fig[1,1])
for _ in 1:4
    traj_sds, = trajectory(sds, 400, [0.5, 0]; Δt = 0.1)
    lines!(ax, traj_sds; color = (:green, 0.25))
    scatter!(ax, traj_sds[1]; color = :green)
end
fig

# We can see that for the stochastic system, even though all trajectories start
# in the basin of the right attractor, over time they will deviate and visit
# the alternative attractor (and then come back and forth forever).

# This is a well known result, for white noise, called "XXX" TODO: .
# And indeed one can study it with CriticalTransitions.jl as we show below!

# ## RandomSystem: example applications

# TODO: probable escale path something something.

# ## RateSystem: example applications

# TODO: Remains to be done.

#