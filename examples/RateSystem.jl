
# # Rate System Example

using CriticalTransitions
using DifferentialEquations
using ModelingToolkit
using CairoMakie
using CairoMakie.Makie.MathTeXEngine: get_font
font = (;
    regular=get_font(:regular),
    bold=get_font(:bold),
    italic=get_font(:italic),
    bold_italic=get_font(:bolditalic),
);

# Let us explore an example of the RateSystem setting of [CriticalTransitions.jl](https://github.com/JuliaDynamics/CriticalTransitions.jl).

# ## Prototypical model for R-tipping, with critical rate r = 4/3

# The following simple one-dimensional model with one attractor is given by the following ordinary differential equations:
# ```math
# \begin{aligned}
#     \dot{x} &= (x+\lambda)^2 - 1
# \end{aligned}
# ```
# The parameter ``\lambda`` shifts the location of the extrema of the drift field.

function f(u, p, t) # out-of-place
    x = u[1]
    λ = p[1]
    dx = (x+λ)^2 - 1
    return SVector{1}(dx)
end
lambda = 0.0
p = [lambda]
x0 = [-1.0]
auto_sys = CoupledODEs(f, x0, p)

## Non-autonomous case

# Now, we define the time-dependent parameter function ``\lambda(t)`` to make the system non-autonomous and investigate the system's behaviour under parameter rampings.

# ```@example RateSystem
# function λ(p,t)
#     λ_max = p[1]
#     lambda = (λ_max/2)*(tanh(λ_max*t/2)+1)
#     return SVector{1}(lambda)
# end
# ```

function λ(p, t)
    λ_max = p[1]
    lambda = (λ_max/2)*(tanh(λ_max*t/2)+1)
    return SVector{1}(lambda)
end

# We define the following parameters
# ```@example RateSystem
# λ_max = 3.
# p_lambda = [λ_max] # parameter of the function lambda
# r = 4/3-0.02 # r just below critical rate
# t_start = -Inf # start time of non-autonomous part
# t_end = Inf    # end time of non-autonomous part
# # And the initial time of the system
# t0 = -10.
# ```

λ_max = 3.0
p_lambda = [λ_max] # parameter of the function lambda
r = 4/3-0.02 # r just below critical rate
t_start = -Inf # start time of non-autonomous part
t_end = Inf    # end time of non-autonomous part
# And the initial time of the system
t0 = -10.0

# We define the RateProtocol

# ```@example RateSystem
# rp = RateProtocol(λ,p_lambda,r,t_start,t_end)
# ```
rp = CriticalTransitions.RateProtocol(λ, p_lambda, r, t_start, t_end)

# We plot the two trajectories

# ```@example RateSystem
# fig = Figure(); axs = Axis(fig[1,1])
# lines!(axs,t0.+auto_traj[2],auto_traj[1][:,1],linewidth=2,label=L"\text{Autonomous system}")
# lines!(axs,nonauto_traj[2],nonauto_traj[1][:,1],linewidth=2,label=L"\text{Nonautonomous system}")
# axislegend(axs,position=:rc,labelsize=10)
# fig
# ```
fig = Figure();
axs = Axis(fig[1, 1])
lines!(
    axs,
    t0 .+ auto_traj[2],
    auto_traj[1][:, 1];
    linewidth=2,
    label=L"\text{Autonomous system}",
)
lines!(
    axs,
    nonauto_traj[2],
    nonauto_traj[1][:, 1];
    linewidth=2,
    label=L"\text{Nonautonomous system}",
)
axislegend(axs; position=:rc, labelsize=10)
fig
