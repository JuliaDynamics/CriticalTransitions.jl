# # [Transition Path Theory for the double well](@id TPT_example)

# In this example, we explore the application of Transition Path Theory (TPT) to a double well system. We will compute various quantities of interest in TPT, such as the Hamiltonian, committor functions, reactive currents, and reaction rates. These computations will be performed on a triangular mesh in the phase space, providing insights into the system's dynamics and transition paths between different states.

# As for now, TPT functionality is considered experimental, and the API may change in the future.

using CriticalTransitions

using CairoMakie
using OrdinaryDiffEq, DelaunayTriangulation, Contour

# ## System

# We consider a simple model for a particle in a double-well potential, subject to dissipation and diffusion. The equation of motion under additive Gaussian noise is given by

# ```math
# \dot{x} = p, \\
# \dot{p} = -\gamma p - \nabla U + \sqrt{\frac{2\gamma}{\beta}} \dot{W},
# ```

# with the potential energy $U(x) = \frac{1}{4}x^4 - \frac{1}{2}x^2$ and the kinetic energy $K(p) = p^2/2$. The parameters $\gamma$ and $\beta=1/k_b T$ control the strength of the dissipation and noise, respectively. $W(t)$ is a Wiener process, and the noise term is scaled by $\sqrt{2\gamma/\beta}$ to ensure the correct temperature scaling for Langevin dynamics defined by the Hamiltonian $H$.

# In CriticalTransitions.jl, a stochastic dynamical system of this form can be constructed as a [`LangevinSystem`](@ref), which takes five input arguments as exemplified below: the Hamiltonian, the divergence-free part of the drift, the system's kinetic energy, the damping coefficient, and the inverse temperature (noise intensity).

beta = 20.0
gamma = 0.5

function Hamiltonian(x, y)
    return 0.5 .* y .^ 2 .+ 0.25 .* x .^ 4 .- 0.5 .* x .^ 2
end

function KE(x)
    return 0.5 .* (x[:, 2] .^ 2)
end

function divfree(x, y)
    f1 = y
    f2 = .-x .^ 3 .+ x
    return f1, f2
end

function double_well(x, y)
    f1 = y
    f2 = .-gamma .* y .- x .^ 3 .+ x
    return f1, f2
end

using CriticalTransitions: LangevinSystem
langevin_sys = LangevinSystem(Hamiltonian, divfree, KE, gamma, beta)

# ## Phase space mesh

# We can easily evaluate and visualize the Hamiltonian on an equally spaced grid in phase space.

nx, ny = 41, 41
nxy = nx * ny
xmin, xmax = -2.0, 2.0
ymin, ymax = -2.0, 2.0

x1 = range(xmin, xmax; length=nx)
y1 = range(ymin, ymax; length=ny)

x_grid = [xx for yy in y1, xx in x1]
y_grid = [yy for yy in y1, xx in x1]

drift1, drift2 = double_well(x_grid, y_grid)
dnorm = sqrt.(drift1 .^ 2 .+ drift2 .^ 2 .+ 1e-12)
H_grid = Hamiltonian(x_grid, y_grid)

fig = CairoMakie.contour(x1, y1, H_grid'; colormap=:viridis, levels=-1:0.4:2, linewidth=2)
v(x::Point2) = Point2f(double_well(x[1], x[2])...)
streamplot!(
    v,
    -2 .. 2,
    -2 .. 2;
    linewidth=0.5,
    colormap=[:black, :black],
    gridsize=(40, 40),
    arrow_size=8,
)
fig

# The double well system is autonomous and respects detailed balance. In this case, the maximum likelihood path (MLP) follows parallel to a flowline of the drift field and can be computed via the string method. If the saddle point is known, one can easily compute the MLP by solving for the (reverse) flow/drift from the saddle point to each of the potential minima. The MLP from (-1,0) to (1,0) gives:

using OrdinaryDiffEq

function reverse_drift!(du, u, p, t)
    du[1] = -u[2]
    return du[2] = -u[2] + u[1] * (u[1]^2 - 1)
end

function drift!(du, u, p, t)
    du[1] = u[2]
    return du[2] = -u[2] - u[1] * (u[1]^2 - 1)
end

prob0 = ODEProblem(reverse_drift!, [-0.001, 0.0], (0.0, 100.0))
sol0 = solve(prob0, Tsit5(); abstol=1e-12, reltol=1e-12)

prob1 = ODEProblem(drift!, [0.001, 0.0], (0.0, 100.0))
sol1 = solve(prob1, Tsit5(); abstol=1e-12, reltol=1e-12)

fig = streamplot(
    v,
    -2 .. 2,
    -2 .. 2;
    linewidth=0.5,
    colormap=[:gray, :gray],
    gridsize=(40, 40),
    arrow_size=8,
)
y = sol0
lines!(y[1, :], y[2, :]; linewidth=2, color=:black)
y = sol1
lines!(y[1, :], y[2, :]; linewidth=2, color=:black)
fig

# Close to the local minima $(-1.0, 0.0)$ and $(1.0, 0.0)$ of the potential landscape, the system under the drift will dissipate to the corresponding attractor. TPT investigates the "reaction" (the name originates from studies of chemical reactions) between two sets in phase space A and B; here we define the two sets to be an ellipse around these minima:

using CriticalTransitions: get_ellipse
point_a = (-1.0, 0.0)
point_b = (1.0, 0.0)
radii = (0.3, 0.4)
density = 0.04

Na = round(Int, π * sum(radii) / density) # the number of points on the A-circle
Nb = Na

ptsA = get_ellipse(point_a, radii, Na)
ptsB = get_ellipse(point_b, radii, Na);

# We also compute an outer boundary of the phase space defined by the maximum value of the Hamiltonian: `H_bound=0.5`. For this, we use the `Contour.jl` package to compute the contour at the level `H_bound`. Just as the ellipse around the attractors, we also re-parameterize the boundary to have a uniform grid spacing.

import Contour as CTR
H_bound = 0.5
cont = CTR.contour(x1, y1, H_grid, H_bound)
yc, xc = coordinates(CTR.lines(cont)[1])
p_outer = [xc yc]

using CriticalTransitions: reparameterization
pts_outer = reparameterization(p_outer, density);
Nouter = size(pts_outer, 1)
Nfix = Na + Nb + Nouter

fig = scatter(ptsA[:, 1], ptsA[:, 2]; label="A points")
scatter!(ptsB[:, 1], ptsB[:, 2]; label="B points")
scatter!(pts_outer[:, 1], pts_outer[:, 2]; label="Outer points")
fig

# We would like to compute the committor, the reactive current, and the reaction rate for the double well with additive Gaussian noise. We compute these quantities on a triangular mesh between the previously computed boundaries.

using CriticalTransitions: distmesh2D, dellipse, ddiff, dunion, huniform
box = [xmin, xmax, ymin, ymax]
pfix = zeros(Nfix, 2)
pfix[1:Na, :] .= ptsA
pfix[(Na + 1):(Na + Nb), :] .= ptsB
pfix[(Na + Nb + 1):Nfix, :] .= pts_outer

function dfunc(p)
    d0 = Hamiltonian(p[:, 1], p[:, 2])
    dA = dellipse(p, point_a, radii)
    dB = dellipse(p, point_b, radii)
    d = ddiff(d0 .- H_bound, dunion(dA, dB))
    return d
end

mesh = distmesh2D(dfunc, huniform, density, box, pfix)

pts, tri = mesh.pts, mesh.tri
fig = Figure()
ax = Axis(fig[1, 1])
for i in 1:size(tri, 1)
    lines!(
        ax,
        [pts[tri[i, j], 1] for j in [1, 2, 3, 1]],
        [pts[tri[i, j], 2] for j in [1, 2, 3, 1]];
        color=:black,
        linewidth=0.1,
    )
end
fig

## Committor functions

# The committor is a scalar function that measures the probability that a system, starting at a given point in phase space, reaches one designated region before another. Formally, for two disjoint sets A and B, the forward committor $q_+(x, p)$ from A to B gives the likelihood that a trajectory initiated at x will reach B before A under the system’s dynamics. The committor boundary-value problem for a Langevin system is given by:

# ```math
#  p \frac{\mathrm{d}q}{\mathrm{d}x} - U'(x) \frac{\mathrm{d}q}{\mathrm{d}p} + \gamma [-p \frac{\mathrm{d}q}{\mathrm{d}p} + \beta^{-1} \frac{\mathrm{d}^2 q}{\mathrm{d}p^2}] = 0,
# ```

# for $(x,p) \in (A\cup B)^c$, with boundary conditions $q(\partial A) = 0$, $q(\partial B) = 1$, and $\nabla \nabla q = 0$ on the outer boundary ${(x,p) : H(x,p) = \mathrm{H_bound}}$. The homogeneous Neumann boundary condition $\nabla \nabla q = 0$ means that the trajectory reflects from the outer boundary whenever it reaches it. We can compute the committor function for the system using the `committor` function.

using CriticalTransitions: committor, find_boundary
_, Aind = find_boundary(mesh.pts, point_a, radii, density)
_, Bind = find_boundary(mesh.pts, point_b, radii, density)

q = committor(langevin_sys, mesh, Aind, Bind)

@show extrema(q)

tricontourf(Triangulation(mesh.pts', mesh.tri'), q)

# We can also compute the backward committor $q_{-}(x, p)$ from A to B, which is the probability that a trajectory initiated at x will reach A before B under the system’s dynamics. Hence, we must reverse the drift function in the Langevin system and swap the boundaries A and B in the committor function

function divfree_reverse(x, y)
    f1, f2 = divfree(x, y)
    return -f1, -f2
end

langevin_sys_reverse = LangevinSystem(Hamiltonian, divfree_reverse, KE, gamma, beta)

qminus = committor(langevin_sys_reverse, mesh, Bind, Aind)

@show extrema(qminus)

tricontourf(Triangulation(mesh.pts', mesh.tri'), qminus)

# For non-equilibrium processes, such as critical transitions in the double-well, we have $q_{-}\neq 1-q_+$. In particular, for Langevin dynamics of the form above, time reversal involves a momentum flip such that $q_{-}(x, p)= 1-q_+(x, -p)$.

# ## Probability density of reactive trajectories

# In general, we are interested in reactive trajectories that start in A and end in B without going back to A. The probability density of finding a reactive trajectory at a point in phase space is given by:

# ```math
# \rho_R(x, p) = \rho(x, p) q₊(x, p) q₋(x, p),
# ```

# where $\rho(x, p)$ is the probability density of finding a trajectory at $(x,p)$, $\rho(x, p)$ is also called the invariant probability density (or invariant measure) of the system. For an overdamped Langevin system the invariant probability density reads:

# ```math
# \rho(x, p) \approx exp(-\beta H(x,p))/Z
# ```

# with $Z=\int exp(-\beta H(x,p)) \mathrm{d}x \mathrm{d}p$ the normalization. We can compute the integrated invariant probability density `Z` for the mesh using the `invariant_pdf` function.

function dfuncA(p)
    return dellipse(p, point_a, radii)
end

function dfuncB(p)
    return dellipse(p, point_b, radii)
end

xa, ya = point_a
xb, yb = point_b
rx, ry = radii
bboxA = [xa - rx, xa + rx, ya - ry, ya + ry]
Amesh = distmesh2D(dfuncA, huniform, density, bboxA, ptsA)
bboxB = [xb - rx, xb + rx, yb - ry, yb + ry]
Bmesh = distmesh2D(dfuncB, huniform, density, bboxB, ptsB)

using CriticalTransitions: invariant_pdf
Z = invariant_pdf(langevin_sys, mesh, Amesh, Bmesh)

@show Z

# Hence, the probability density of a reactive trajectory is given by:

mu = exp.(-beta * Hamiltonian(pts[:, 1], pts[:, 2])) / Z
muAB = mu .* q .* qminus

tricontourf(Triangulation(mesh.pts', mesh.tri'), muAB)

# The current of reactive trajectories is given by:

# ```math
# J_R = \frac{e^{-\beta H} q_{+} q_{-}}{Z}\binom{p}{-\nabla U}+k_B T \gamma Z^{-1} e^{-\beta H}\binom{0}{q_{-} \frac{\partial q_{+}}{\partial p}-q_{+} \frac{\partial q_{-}}{\partial p}}
# ```

# and the transition rate:

# ```math
# v_R=k_B T \gamma Z_H^{-1} \int \sum_{i=1}^d m_i\left(\frac{\partial q_{+}}{\partial p_i}\right)^2 e^{-\beta H({x}, p)} d {x} d {p}
# ```

# These can be computed using the `reactive_current` function:

using CriticalTransitions: reactive_current
Rcurrent, Rrate = reactive_current(langevin_sys, mesh, q, qminus, Z)
@show Rrate

# Plotting the current norm reveals that the current is the strongest around the saddle point.

ARcurrent = vec(sqrt.(sum(Rcurrent .^ 2; dims=2)))
ARCmax = maximum(ARcurrent)

tricontourf(Triangulation(mesh.pts', mesh.tri'), ARcurrent)

# The transition current has a direction from A to B.

c = ARcurrent ./ ARCmax
arrows2d(
    pts[:, 1],
    pts[:, 2],
    Rcurrent[:, 1] ./ ARCmax,
    Rcurrent[:, 2] ./ ARCmax;
    color=c,
    lengthscale=0.1,
)

#

using CriticalTransitions: probability_reactive
prob_reactive = probability_reactive(langevin_sys, mesh, q, qminus, Z)
print(
    "Probability that a trajectory is reactive at a randomly picked time: ", prob_reactive
)

#
using CriticalTransitions: probability_last_A
prob_lastA = probability_last_A(langevin_sys, mesh, Amesh, qminus, Z)
print("Probability that a trajectory last visited A: ", prob_lastA)

# Authored by O. Ameye and R. Börner
