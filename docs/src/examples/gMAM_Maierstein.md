
# Minimal action path for the Maier-Stein model.

Here we compute the maximum likelihood transition path (instanton) for the Maier-Stein model.

```@example GMAM
using CriticalTransitions

using CairoMakie
using CairoMakie.Makie.MathTeXEngine: get_font
font = (;
    regular=get_font(:regular), bold=get_font(:bold),
    italic=get_font(:italic), bold_italic=get_font(:bolditalic)
)
```

Let us explore the features of [CriticalTransitions.jl](https://github.com/JuliaDynamics/CriticalTransitions.jl) with Maier-Stein model.

## Maier-stein model

The [Maier-Stein model](https://link.springer.com/article/10.1007/BF02183736) (J. Stat. Phys 83, 3–4 (1996)) is commonly used in the field of nonlinear dynamics for benchmarking Large Deviation Theory (LDT) techniques, e.g., stoachastic transitions between different stable states. It is a simple model that describes the dynamics of a system with two degrees of freedom ``u`` and ``v``, and is given by the following set of ordinary differential equations:
```math
\begin{aligned}
    \dot{u} &= u-u^3 - \beta*u*v^2\\
    \dot{v} &= -\alpha (1+u^2)*v
\end{aligned}
```
The parameter ``\alpha>0`` controls the strength of the drift field and ``\beta>0`` represents the softening of that drift field.

```@example GMAM
function meier_stein!(du, u, p, t) # in-place
    x, y = u
    du[1] = x - x^3 - 10 * x * y^2
    du[2] = -(1 + x^2) * y
end
function meier_stein(u, p, t) # out-of-place
    x, y = u
    dx = x - x^3 - 10 * x * y^2
    dy = -(1 + x^2) * y
    SA[dx, dy]
end
σ = 0.25
sys = StochSystem(meier_stein, [], zeros(2), σ, idfunc, nothing, I(2), "WhiteGauss")
```

A good reference to read about the large deviations methods is [this](https://homepages.warwick.ac.uk/staff/T.Grafke/simplified-geometric-minimum-action-method-for-the-computation-of-instantons.html) or [this]( https://homepages.warwick.ac.uk/staff/T.Grafke/simplified-geometric-minimum-action-method-for-the-computation-of-instantons.html) blog post by Tobias Grafke.

## Attractors

We start by investigating the deterministic dynamics of the Maier-Stein model.

The function `fixed points` return the attractors, their eigenvalues and stability within the state space volume defined by `bmin` and `bmax`.

```@example GMAM
u_min = -1.1;
u_max = 1.1;
v_min = -0.4;
v_max = 0.4;
bmin = [u_min, v_min];
bmax = [u_max, v_max];
fp, eig, stab = fixedpoints(sys, bmin, bmax)
stable_fp = fp[stab]
```

```@example GMAM
res = 100
u_range = range(u_min, u_max, length=res)
v_range = range(v_min, v_max, length=res)

du(u, v) = u - u^3 - 10 * u * v^2
dv(u, v) = -(1 + u^2) * v
odeSol(u, v) = Point2f(du(u, v), dv(u, v))

z = [norm([du(x, y), dv(x, y)]) for x in u_range, y in v_range]
zmin, zmax = minimum(z), maximum(z)

fig = Figure(size=(600, 400), fontsize=13)
ax = Axis(fig[1, 1], xlabel="u", ylabel="v", aspect=1.4,
    xgridcolor=:transparent, ygridcolor=:transparent,
    ylabelrotation=0)

hm = heatmap!(ax, u_range, v_range, z, colormap=:Blues, colorrange=(zmin, zmax))
Colorbar(fig[1, 2], hm; label="", width=15, ticksize=15, tickalign=1)
streamplot!(ax, odeSol, (u_min, u_max), (v_min, v_max);
    gridsize=(20, 20), arrow_size=10, stepsize=0.01,
    colormap=[:black, :black]
)
colgap!(fig.layout, 7)
limits!(u_min, u_max, v_min, v_max)
fig

[scatter!(ax, Point(fp[i]), color=stab[i] > 0 ? :red : :dodgerblue,
    markersize=10) for i in eachindex(fp)]
fig
```

We can simulate a stochastic trajectory using the function `simulate`.

```@example GMAM
sol = simulate(sys, [1.0, 0.0], tmax=1000)

fig = Figure(size=(1000, 400), fontsize=13)
ax1 = Axis(fig[1, 1], xlabel="t", ylabel="u", aspect=1.2,
    xgridcolor=:transparent, ygridcolor=:transparent,
    ylabelrotation=0)
ax2 = Axis(fig[1, 2], xlabel="u", ylabel="v", aspect=1.2,
    xgridcolor=:transparent, ygridcolor=:transparent,
    ylabelrotation=0)

lines!(ax1, sol.t, first.(sol.u); linewidth=2, color=:black)

hm = heatmap!(ax2, u_range, v_range, z, colormap=:Blues, colorrange=(zmin, zmax))
Colorbar(fig[1, 3], hm; label="", width=15, ticksize=15, tickalign=1)
streamplot!(ax2, odeSol, (u_min, u_max), (v_min, v_max);
    gridsize=(20, 20), arrow_size=10, stepsize=0.01,
    colormap=[:white, :white]
)
colgap!(fig.layout, 7)
limits!(u_min, u_max, v_min, v_max)
fig

[scatter!(ax2, Point(fp[i]), color=stab[i] > 0 ? :red : :dodgerblue,
    markersize=10) for i in eachindex(fp)]

lines!(ax2, reduce(hcat, sol.u), linewidth=1, color=(:black, 0.2))
fig
```

## Basins of attraction

Basins of attraction are the regions in the state space that lead to a particular attractor. We can find the basins of attraction using the function `basins`.

@example GMAM
ba = basins(sys, [0.0, 0], [0.0, 1], [1.0, 0], intervals_to_box([-2, -2], [2, 2]), bstep=[0.01, 0.01], ϵ_mapper=0.001, Ttr=100)
Ur, Vr, atr, M = ba
heatmap(Ur, Vr, M)
```

The basin boundaries can be quickly extracted using the function `basin_boundaries`.

@example GMAM
bb = basinboundary(ba)
```

## Edge tracking

The edge tracking algorithm is a simple numerical method to find the edge state or (possibly chaotic) saddle on the boundary between two basins of attraction. It is first introduced by [Battelino et al. (1988)](https://doi.org/10.1016/0167-2789(88)90057-7) and further described by [Skufca et al. (2006)](https://doi.org/10.1103/PhysRevLett.96.174101).

@example GMAM
edge = edgetracking(sys, [0.0, 0.0], [-1.0, -0.2], [stable_fp[1], stable_fp[2]])
```

## Transitions

We can quickly find a path which computes a transition from one attractor to another using the function `transition.

@example GMAM
path, time, succes = transition(sys, fp[stab][1], fp[stab][2])
```

@example GMAM
fig = Figure(size=(600, 400), fontsize=13)
ax = Axis(fig[1, 1], xlabel="u", ylabel="v", aspect=1.4,
    xgridcolor=:transparent, ygridcolor=:transparent,
    ylabelrotation=0)

hm = heatmap!(ax, u_range, v_range, z, colormap=:Blues, colorrange=(zmin, zmax))
Colorbar(fig[1, 2], hm; label="", width=15, ticksize=15, tickalign=1)
streamplot!(ax, odeSol, (u_min, u_max), (v_min, v_max);
    gridsize=(20, 20), arrow_size=10, stepsize=0.01,
    colormap=[:white, :white]
)
colgap!(fig.layout, 7)
limits!(u_min, u_max, v_min, v_max)
fig

[scatter!(ax, Point(fp[i]), color=stab[i] > 0 ? :red : :dodgerblue,
    markersize=10) for i in eachindex(fp)]
fig

lines!(ax, path.u, lw=2, color=:black)
fig
```

If we want to compute many: `transitions` is the function to use.

@example GMAM
tt = transitions(sys, fp[stab][1], fp[stab][2], 1, rad_i=0.1, rad_f=0.1, tmax=1e3);
```

@example GMAM
fig = Figure(size=(600, 400), fontsize=13)
ax = Axis(fig[1, 1], xlabel="u", ylabel="v", aspect=1.4,
    xgridcolor=:transparent, ygridcolor=:transparent,
    ylabelrotation=0)

hm = heatmap!(ax, u_range, v_range, z, colormap=:Blues, colorrange=(zmin, zmax))
Colorbar(fig[1, 2], hm; label="", width=15, ticksize=15, tickalign=1)
streamplot!(ax, odeSol, (u_min, u_max), (v_min, v_max);
    gridsize=(20, 20), arrow_size=10, stepsize=0.01,
    colormap=[:black, :black]
)
colgap!(fig.layout, 7)
limits!(u_min, u_max, v_min, v_max)
fig

[scatter!(ax, Point(fp[i]), color=stab[i] > 0 ? :red : :dodgerblue,
    markersize=10) for i in eachindex(fp)]

for i in 1:3
    lines!(ax, tt[1][i].u)
end
fig
```

@example GMAM
# R, L = fp[stab]
# initial = reduce(hcat, path.u)
# lv = langevinmcmc(sys, initial; tmax =0.1)
```

## Large deviation theory

In the context of nonlinear dynamics, Large Deviation Theory provides tools to quantify the probability of rare events that deviate significantly from the system's typical behavior. These rare events might be extreme values of a system's output, sudden transitions between different states, or other phenomena that occur with very low probability but can have significant implications for the system's overall behavior.

Large deviation theory applies principles from probability theory and statistical mechanics to develop a rigorous mathematical description of these rare events. It uses the concept of a rate function, which measures the exponential decay rate of the probability of large deviations from the mean or typical behavior. This rate function plays a crucial role in quantifying the likelihood of rare events and understanding their impact on the system.

For example, in a system exhibiting chaotic behavior, LDT can help quantify the probability of sudden large shifts in the system's trajectory. Similarly, in a system with multiple stable states, it can provide insight into the likelihood and pathways of transitions between these states under fluctuations. In the context of the Minimum Action Method (MAM) and the Geometric Minimum Action Method (gMAM), Large Deviation Theory is used to handle the large deviations action functional on the space of curves. This is a key part of how these methods analyze dynamical systems.

The Maier-Stein model is a typical benchmark to test such LDT techniques. Let us try to reproduce the following figure from [Tobias Grafke's blog post](https://homepages.warwick.ac.uk/staff/T.Grafke/rogue-waves-and-large-deviations.html):

![image-2.png](attachment:image-2.png)

Let us first make an initial path:

@example GMAM
xx = range(-1.0, 1.0, length=100)
yy = 0.3 .* (-xx .^ 2 .+ 1)
init = Matrix([xx yy]')
```

`min_action_method` runs the Minimum Action Method (MAM) to find the minimum action path (instanton) between an initial state x_i and final state x_f. This algorithm uses the minimizers of the [Optim.jl](https://julianlsolvers.github.io/Optim.jl/stable/#) package to minimize the [Freidlin-Wentzell action functional](https://en.wikipedia.org/wiki/Freidlin%E2%80%93Wentzell_theorem) (see `fw_action`) or [Onsager-Machlup action functional](https://en.wikipedia.org/wiki/Onsager%E2%80%93Machlup_function) (see `om_action`) for the given StochSystem `sys`.

@example GMAM
mm = min_action_method(sys, init, 1, maxiter=200, save_info=false, verbose=false)
```

@example GMAM
fig = Figure(size=(600, 400), fontsize=13)
ax = Axis(fig[1, 1], xlabel="u", ylabel="v", aspect=1.4,
    xgridcolor=:transparent, ygridcolor=:transparent,
    ylabelrotation=0)

hm = heatmap!(ax, u_range, v_range, z, colormap=:Blues, colorrange=(zmin, zmax))
Colorbar(fig[1, 2], hm; label="", width=15, ticksize=15, tickalign=1)
streamplot!(ax, odeSol, (u_min, u_max), (v_min, v_max);
    gridsize=(20, 20), arrow_size=10, stepsize=0.01,
    colormap=[:black, :black]
)
colgap!(fig.layout, 7)
limits!(u_min, u_max, v_min, v_max)
fig

[scatter!(ax, Point(fp[i]), color=stab[i] > 0 ? :red : :dodgerblue,
    markersize=10) for i in eachindex(fp)]

lines!(ax, init, linewidth=3, color=:black, linestyle=:dash)
lines!(ax, mm, linewidth=3, color=:orange)
fig
```

TODO: Why does MAM get the wrong results?

`geometric_min_action_method` computes the minimizer of the Freidlin-Wentzell action using the geometric minimum action method (gMAM). The Minimum Action Method (MAM) is a more traditional approach, while the Geometric Minimum Action Method (gMAM) is a blend of the original MAM and the [string method](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.66.052301).

@example GMAM
gm = geometric_min_action_method(sys, init, converge=1e-7, maxiter=1000)
```

@example GMAM
fig = Figure(size=(600, 400), fontsize=13)
ax = Axis(fig[1, 1], xlabel="u", ylabel="v", aspect=1.4,
    xgridcolor=:transparent, ygridcolor=:transparent,
    ylabelrotation=0)

hm = heatmap!(ax, u_range, v_range, z, colormap=:Blues, colorrange=(zmin, zmax))
Colorbar(fig[1, 2], hm; label="", width=15, ticksize=15, tickalign=1)
streamplot!(ax, odeSol, (u_min, u_max), (v_min, v_max);
    gridsize=(20, 20), arrow_size=10, stepsize=0.01,
    colormap=[:black, :black]
)
colgap!(fig.layout, 7)
limits!(u_min, u_max, v_min, v_max)
fig

[scatter!(ax, Point(fp[i]), color=stab[i] > 0 ? :red : :dodgerblue,
    markersize=10) for i in eachindex(fp)]

# len = length(gm[1])
# for i in 2:(len-20)
#     i & 10 == 0 && lines!(ax, gm[1][i], linewidth = 3, color = (:black, i/(3*len)))
# end
lines!(ax, gm[1][1], linewidth=3, color=:black, linestyle=:dash)
lines!(ax, gm[1][end], linewidth=3, color=:orange)
# lines!(ax, mm, linewidth = 3, color = :orange)
fig
```

---
Author: Orjan Ameye (orjan.ameye@hotmail.com)
Date: 13 Feb 2024