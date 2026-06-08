# # Quasipotential of the Maier-Stein model

using CriticalTransitions
using CairoMakie

# The Freidlin-Wentzell quasipotential ``U_A(x)`` measures the minimal action of the path connecting an attractor ``x_A`` to a point ``x`` in the small-noise limit. It controls the exponential rate of escape from ``x_A`` and, for non-gradient drift, replaces the usual potential as the relevant Lyapunov function for the deterministic flow plus noise.

# The Ordered Line Integral Method ([Dahiya and Cameron 2018](https://doi.org/10.1007/s10915-017-0590-9)) computes ``U_A(x)`` on a Cartesian grid in a single Dijkstra-style sweep. We apply it to the non-gradient Maier-Stein drift
# ```math
# \begin{aligned}
#     \dot{u} &= u - u^3 - \beta\, u v^2 \\
#     \dot{v} &= -(1 + u^2) v
# \end{aligned}
# ```
# the standard nongradient benchmark in the large-deviations literature. For every ``\beta > 0`` the system has two stable fixed points at ``(\pm 1, 0)`` and a saddle at the origin; the saddle value of ``U_A`` is the Freidlin-Wentzell barrier between the wells.

maierstein(β) = (x, p, t) -> SVector(x[1] - x[1]^3 - β * x[1] * x[2]^2, -(1 + x[1]^2) * x[2])
grid = CartesianGrid((-1.5, 1.5, 121), (-1.0, 1.0, 121))

# ## An exact reference value
#
# On the invariant ``u``-axis (``v = 0``) the drift collapses to the gradient system ``\dot u = u - u^3`` with potential ``V(u) = -u^2/2 + u^4/4``. As long as the minimum-action path stays on the axis, the quasipotential is therefore known in closed form,
# ```math
# U_A(u, 0) = 2\,\bigl(V(u) - V(-1)\bigr) = \tfrac{1}{2}(1 - u^2)^2 ,
# ```
# with a barrier ``U_A(0,0) = 1/2`` that is *independent of* ``\beta``. This holds while ``\beta`` stays below the focusing threshold ``\beta_c = 4`` ([Maier and Stein 1993](https://doi.org/10.1103/PhysRevE.48.931)); above it the optimal escape path bows off the axis and ``U_A`` is no longer differentiable.
#
# We validate the solver against this exact value at ``\beta = 2``:

saddle = CartesianIndex(argmin(abs.(grid.centers[1])), argmin(abs.(grid.centers[2])))

sys = CoupledSDEs(maierstein(2.0), [-1.0, 0.0]; noise_strength = 0.3)
qp = quasipotential(sys, grid, [-1.0, 0.0])
qp.U[saddle]   # ≈ 0.5

# Overlaying the computed on-axis profile ``U_A(u, 0)`` (the ``v = 0`` row) on the closed form confirms the match on the source side ``u \in [-1, 0]``. Past the saddle the system descends deterministically for free, so ``U_A`` stays flat at ``1/2`` and the formula (dashed) no longer applies:

us = grid.centers[1]
jaxis = saddle[2]   # grid row closest to v = 0

fig, ax, _ = lines(us, qp.U[:, jaxis]; label = "OLIM  U_A(u, 0)", axis = (; xlabel = "u", ylabel = "U_A(u, 0)"))
uleft = range(-1, 0; length = 100)
lines!(ax, uleft, @. (1 - uleft^2)^2 / 2; linestyle = :dash, label = "½(1 - u²)²")
axislegend(ax; position = :ct)
fig

# ## Off-axis regime ``\beta = 5``
#
# For ``\beta = 5 > \beta_c`` the minimum-action path leaves the axis and the barrier drops below ``1/2``:

sys = CoupledSDEs(maierstein(5.0), [-1.0, 0.0]; noise_strength = 0.3)
qp = quasipotential(sys, grid, [-1.0, 0.0])
qp.U[saddle]   # ≈ 0.483

# Contour plot of `qp.U`, with both attractors marked in white and the saddle in red.

fig, ax, _ = contourf(
    qp.grid.centers[1], qp.grid.centers[2], qp.U;
    levels = 30, colormap = :viridis,
)
scatter!(
    ax, [-1.0, 1.0, 0.0], [0.0, 0.0, 0.0];
    color = [:white, :white, :red], markersize = 14,
)
fig
