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
# with ``\beta = 5``. The system has two stable fixed points at ``(\pm 1, 0)`` and a saddle at the origin; the saddle value of ``U_A`` is the Freidlin-Wentzell barrier between the wells.

f(x, p, t) = SVector(x[1] - x[1]^3 - 5 * x[1] * x[2]^2, -(1 + x[1]^2) * x[2])
sys = CoupledSDEs(f, [-1.0, 0.0]; noise_strength = 0.3)

# We pick a grid that covers both attractors and call `quasipotential` with the left attractor as the source.

grid = CartesianGrid((-1.5, 1.5, 121), (-1.0, 1.0, 81))
qp = quasipotential(sys, grid, [-1.0, 0.0])

# The saddle at ``(0, 0)`` sits at cell `CartesianIndex(61, 41)`. For this system the analytic barrier is known to be ``1/3``.

qp.U[CartesianIndex(61, 41)]

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
