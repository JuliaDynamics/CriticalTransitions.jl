# Transition Path Theory for the undriven Duffing oscillator

In this example, we explore the application of Transition Path Theory to the undriven Duffing oscillator. We will compute various quantities of interest in Transition Path Theory (TPT), such as the Hamiltonian, committor functions, reactive currents, and reaction rates. These computations will be performed on a triangular mesh in the phase space, providing insights into the system's dynamics and transition paths between different states.

```@example TPT
using CriticalTransitions

using CairoMakie
using OrdinaryDiffEq, DelaunayTriangulation, Contour
```

## System

The Duffing oscillator is a simple model for a nonlinear oscillator with a double-well potential. The equation of motion for the Duffing oscillator under additive Gaussian noise is given by:

```math
\dot{x} = p, \\
\dot{p} = -\gamma p - \nabla U + \sqrt{\frac{2\gamma}{\beta}} \dot{W},
```

with the potential energy $U(x) = \frac{1}{4}x^4 - \frac{1}{2}x^2$ and the kinetic energy $K(p) = p^2/2$. The parameters $\gamma$ and $\beta=1/k_b T$ control the strength of the dissipation and noise, respectively. $W$ is a Wiener process, and the noise term is scaled by $\sqrt{2\gamma/\beta}$ to ensure the correct temperature scaling for a Langevin type system defined by the Hamiltonian $H$.

```@example TPT
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

function duffing(x, y)
    f1 = y
    f2 = .- gamma .* y .- x.^3 .+ x
    return f1, f2
end

langevin_sys = Langevin(Hamiltonian, divfree, KE, gamma, beta)
```

## Phase space mesh

We can easily evaluate and visualize the Hamiltonian and equally spaced grid in phase space.

```@example TPT
nx, ny = 41, 41
nxy = nx * ny
xmin, xmax = -2.0, 2.0
ymin, ymax = -2.0, 2.0

x1 = range(xmin, xmax, length=nx)
y1 = range(ymin, ymax, length=ny)

x_grid = [xx for yy in y1, xx in x1]
y_grid = [yy for yy in y1, xx in x1]

drift1, drift2 = duffing(x_grid, y_grid)
dnorm = sqrt.(drift1.^2 .+ drift2.^2 .+ 1e-12)
Hgrid = Hamiltonian(x_grid, y_grid)
```

```@example TPT
fig = CairoMakie.contour(x1, y1, Hgrid', colormap = :viridis, levels=-1:0.4:2, linewidth = 2)
v(x::Point2) = Point2f(duffing(x[1], x[2])...)
streamplot!(v, -2..2, -2..2, linewidth = 0.5, colormap = [:black, :black], gridsize = (40, 40), arrow_size = 8)
fig
```

The undriven duffing oscillator, is autonomous and respect detailed balance. As such the maximum likelihood path is the path that is parallel to drift and can be computed with the string method. If one know the saddle point, one can easily compute the MLP by solving for the (reverse) flow/drift from the saddle point to the minima. As such, the maximum likelihood transition path from (-1,0) to (1,0) gives:

```@example TPT
# Generate and plot the maximum likelihood transition path from (-1,0) to (1,0)
using OrdinaryDiffEq

function reverse_drift!(du, u, p, t)
    du[1] = -u[2]
    du[2] = -u[2] + u[1]*(u[1]^2 - 1)
end

function drift!(du, u, p, t)
    du[1] =  u[2]
    du[2] = -u[2] - u[1]*(u[1]^2 - 1)
end

prob0 = ODEProblem(reverse_drift!, [-0.001, 0.0], (0.0, 100.0))
sol0 = solve(prob0, Tsit5(); abstol=1e-12, reltol=1e-12)

prob1 = ODEProblem(drift!, [0.001, 0.0], (0.0, 100.0))
sol1 = solve(prob1, Tsit5(); abstol=1e-12, reltol=1e-12)

fig = streamplot(v, -2..2, -2..2, linewidth = 0.5, colormap = [:gray, :gray], gridsize = (40, 40), arrow_size = 8)
y = sol0
lines!(y[1,:],y[2,:],linewidth = 2, color = :black)
y = sol1
lines!(y[1,:],y[2,:],linewidth = 2, color = :black)
fig
```

We have two minima in the potential landscape, such that the system under the drift will dissipate to these corresponding attractors close to $(-1.0, 0.0)$ and $(1.0, 0.0)$. Transition path theory investigates the "reaction" between two sets in phase space A and B, as such we define the two sets to be an ellipse around these minima:

```@example TPT
point_a = (-1.0, 0.0)
point_b = (1.0, 0.0)
radii = (0.3, 0.4)
density = 0.04

Na = round(Int, π * sum(radii) / density) # the number of points on the A-circle
Nb = Na

ptsA = get_ellipse(point_a, radii, Na)
ptsB = get_ellipse(point_b, radii, Na);
```

We also compute an outer boundary of the phase space defined by the maximum value of the Hamiltonian: `Hbdry=0.5`. For this, we use the contour package to compute the contour at the level `Hbdry`. Just as the ellipse around the attractors, we also reparametrize the boundary to have a uniform grid spacing.

```@example TPT
import Contour as CTR
Hbdry = 0.5
cont = CTR.contour(x1, y1, Hgrid, Hbdry)
yc, xc = coordinates(CTR.lines(cont)[1])
p_outer = [xc yc]

pts_outer = reparametrization(p_outer,density);
Nouter = size(pts_outer, 1)
Nfix = Na+Nb+Nouter

fig = scatter(ptsA[:,1], ptsA[:,2], label="A points")
scatter!(ptsB[:,1], ptsB[:,2], label="B points")
scatter!(pts_outer[:,1], pts_outer[:,2], label="Outer points")
fig
```

We would like to compute the committor, the reactive current, and the reaction rate for the Duffing oscillator with additive Gaussian noise. We compute these quantities on a triangular mesh between the before computed boundaries.

```@example TPT
box = [xmin, xmax, ymin, ymax]
pfix = zeros(Nfix, 2)
pfix[1:Na, :] .= ptsA
pfix[Na+1:Na+Nb, :] .= ptsB
pfix[Na+Nb+1:Nfix, :] .= pts_outer

function dfunc(p)
    d0 = Hamiltonian(p[:, 1], p[:, 2])
    dA = dellipse(p, point_a, radii)
    dB = dellipse(p, point_b, radii)
    d = ddiff(d0 .- Hbdry, dunion(dA, dB))
    return d
end

mesh = distmesh2D(dfunc, huniform, density, box, pfix)

pts, tri = mesh.pts, mesh.tri
fig = Figure()
ax = Axis(fig[1, 1])
for i in 1:size(tri,1)
    lines!(ax, [pts[tri[i,j],1] for j in [1,2,3,1]],
          [pts[tri[i,j],2] for j in [1,2,3,1]],
          color=:black, linewidth=0.1)
end
fig
```

## Committor functions

A committor measures the probability that a system, starting at a given point in phase space, will reach one designated region before another. Formally, for two disjoint sets A and B, the forward committor $q_+(x, p)$ from A to B gives the likelihood that a trajectory initiated at x will reach B before A under the system’s dynamics. The committor boundary-value problem for a Langevin system is given by:

```math
 p \frac{\mathrm{d}q}{\mathrm{d}x} - U'(x) \frac{\mathrm{d}q}{\mathrm{d}p} + \gamma [-p \frac{\mathrm{d}q}{\mathrm{d}p} + \beta^{-1} \frac{\mathrm{d}^2 q}{\mathrm{d}p^2}] = 0,
```

for $(x,p) \in (A\cup B)^c$, with boundary conditions $q(\partial A) = 0$, $q(\partial B) = 1$, and $\nabla \nabla q = 0$ on the outer boundary ${(x,p) : H(x,p) = \mathrm{Hbdry}}$. The homogeneous Neumann boundary condition $\nabla \nabla q = 0$ means that the trajectory reflects from the outer boundary whenever it reaches it. We can compute the committor function for the Duffing oscillator using the `committor` function.

```@example TPT
_, Aind = find_boundary(mesh.pts, point_a, radii, density)
_, Bind = find_boundary(mesh.pts, point_b, radii, density)

q = committor(langevin_sys, mesh, Aind, Bind)

@show extrema(q)

tricontourf(Triangulation(mesh.pts', mesh.tri'), q)
```

We can also compute the backward committor $q_{-}(x, p)$ from A to B, which is the probability that a trajectory initiated at x will reach A before B under the system’s dynamics. Hence, we must reverse the drift function in the Langevin system and swap the boundaries A and B in the committor function

```@example TPT
function divfree1(x,y)
    f1,f2 = divfree(x,y)
    return -f1,-f2
end

langevin_sys_reverse = CriticalTransitions.Langevin(Hamiltonian, divfree1, KE, gamma, beta)

qminus = committor(langevin_sys_reverse, mesh, Bind, Aind)

@show extrema(qminus)

tricontourf(Triangulation(mesh.pts', mesh.tri'), qminus)
```

For non-equilibrium processes, such as the transitions in the double-well of the Duffing, we have that the $q_{-}\neq 1-q_+$. In particular, for Langevin systems of the form in the system above time reversal involves a momentum flip such that $q_{-}(x, p)= 1-q_+(x, -p)$.

## Probability Density of Reactive Trajectories

In general, we are interested in reactive trajectories that start in A and ends in B without going back to A. The probability density of finding a reactive trajectory at a point in phase space is given by:

```math
\rho_R(x, p) = \rho(x, p) q₊(x, p) q₋(x, p),
```

where $\rho(x, p)$ is the probability density of finding a trajectory at $(x,p)$, $\rho(x, p)$ is also called the invariant probability density of the system. For an overdamped Langevin system the invariant probability density:

```math
\rho(x, p) \approx exp(-\beta H(x,p))/Z
```

with $Z=\int exp(-\beta H(x,p)) \mathrm{d}x \mathrm{d}p$ the normalization. We can compute the integrated invariant probability density `Z` for the mesh using the `invariant_pdf` function.

```@example TPT
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

Z = invariant_pdf(langevin_sys, mesh, Amesh, Bmesh)

@show Z
```

Hence, the probability density of a reactive trajectory is given by:

```@example TPT
# probability density of reactive trajectories
mu = exp.(-beta * Hamiltonian(pts[:,1], pts[:,2])) / Z
muAB = mu .* q .* qminus

tricontourf(Triangulation(mesh.pts', mesh.tri'), muAB)
```

The current of reactive trajectories is given by:

```math
J_R = \frac{e^{-\beta H} q_{+} q_{-}}{Z}\binom{p}{-\nabla U}+k_B T \gamma Z^{-1} e^{-\beta H}\binom{0}{q_{-} \frac{\partial q_{+}}{\partial p}-q_{+} \frac{\partial q_{-}}{\partial p}}
```

and the transition rate:

```math
v_R=k_B T \gamma Z_H^{-1} \int \sum_{i=1}^d m_i\left(\frac{\partial q_{+}}{\partial p_i}\right)^2 e^{-\beta H({x}, p)} d {x} d {p}
```

These can be computed using the `reactive_current` function.

```@example TPT
Rcurrent, Rrate = reactive_current(langevin_sys, mesh, q, qminus, Z)
@show Rrate
```

Plotting the current norm reveals that the current is the strongest around the saddle point.

```@example TPT
ARcurrent = vec(sqrt.(sum(Rcurrent.^2, dims=2)))
ARCmax = maximum(ARcurrent)

tricontourf(Triangulation(mesh.pts', mesh.tri'), ARcurrent)
```

The transition current has a direction from A to B.

```@example TPT
c =ARcurrent./ARCmax
arrows(pts[:,1], pts[:,2], Rcurrent[:,1]./ARCmax, Rcurrent[:,2]./ARCmax, arrowsize = c*10, lengthscale = 0.1, arrowcolor = c, linecolor = c)
```

```@example TPT
prob_reactive = probability_reactive(langevin_sys, mesh, q, qminus, Z)
print("Probability that a trajectory is reactive at a randomly picked time: ",prob_reactive)
```

```@example TPT
prob_lastA = probability_last_A(langevin_sys, mesh, Amesh, qminus, Z)
print("Probability that a trajectory last visited A: ", prob_lastA)
```
