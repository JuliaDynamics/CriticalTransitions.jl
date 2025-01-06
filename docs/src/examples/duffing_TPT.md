# Transition Path Theory (TPT) for the undriven Duffing oscillator

```@example
using CriticalTransitions

using CairoMakie
using OrdinaryDiffEq, DelaunayTriangulation, Contour
```

```@example
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

function drift(x, y)
    f1 = y
    f2 = .- gamma .* y .- x.^3 .+ x
    return f1, f2
end

langevin_sys = CriticalTransitions.Langevin(Hamiltonian, divfree, KE, gamma, beta)
```

We can easily evaluate the Hamiltonian and equally spaced grid in phase space.

```@example
nx, ny = 41, 41
nxy = nx * ny
xmin, xmax = -2.0, 2.0
ymin, ymax = -2.0, 2.0

x1 = range(xmin, xmax, length=nx)
y1 = range(ymin, ymax, length=ny)

x_grid = [xx for yy in y1, xx in x1]
y_grid = [yy for yy in y1, xx in x1]

drift1, drift2 = drift(x_grid, y_grid)
dnorm = sqrt.(drift1.^2 .+ drift2.^2 .+ 1e-12)
Hgrid = Hamiltonian(x_grid, y_grid)
```

```@example
fig = CairoMakie.contour(x1, y1, Hgrid', colormap = :viridis, levels=-1:0.4:2, linewidth = 2)
v(x::Point2) = Point2f(drift(x[1], x[2])...)
streamplot!(v, -2..2, -2..2, linewidth = 0.5, colormap = [:black, :black], gridsize = (40, 40), arrow_size = 8)
fig
```

We have two minima in the potential landscape, such that the system under the drift will dissipate to these corresponding atractors at close to $(-1.0, 0.0)$ and $(1.0, 0.0)$. We compute two ellepses around these minima:

```@example
point_a = (-1.0, 0.0)
point_b = (1.0, 0.0)
radii = (0.3, 0.4)
density = 0.1

Na = round(Int, π * sum(radii) / density) # the number of points on the A-circle
Nb = Na

ptsA = get_ellipse(point_a, radii, Na)
ptsB = get_ellipse(point_b, radii, Na);
```

We also compute an outer boundary of the phase space defined by the maximum value of the Hamiltonian: `Hbdry=0.5`. For this, we us the contour package to comute the contour at the level `Hbdry`. Just as the ellipse around the attractors, we also reparmetrize the boundary to have a uniform grid spacing.

```@example
Hbdry = 0.5
cont = Contour.contour(x1, y1, Hgrid, Hbdry)
yc, xc = coordinates(Contour.lines(cont)[1])
p_outer = [xc yc]

pts_outer = reparametrization(p_outer,density);
Nouter = size(pts_outer, 1)
Nfix = Na+Nb+Nouter
```

We plot the computed boundaries, for which we compute a triangulation using the distmesh module.

```@example
fig = scatter(ptsA[:,1], ptsA[:,2], label="A points")
scatter!(ptsB[:,1], ptsB[:,2], label="B points")
scatter!(pts_outer[:,1], pts_outer[:,2], label="Outer points")
fig
```

The overdamped undriven duffing oscillator, is autonomous and respect detailed balance. As such the maximum likelihood path is the path that is parallel to drift and can easily be computed with the string method. If one know the sadde point, one can easily compute the MLP by solving for the (reverse) flow/drift from the saddle point to the minima. As such,  the maximum likelihood transition path from (-1,0) to (1,0) gives:

```@example
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

# plot!(legend=:topright, aspect_ratio=:equal)
y = sol0
lines!(y[1,:],y[2,:],linewidth = 2, color = :black)
y = sol1
lines!(y[1,:],y[2,:],linewidth = 2, color = :black)
fig
```

We would like to compute the committor, the reactive current, and the reaction rate for the overdamped Duffing oscillator with additive gaussian noise. We compute this quantities on a triangular mesh between the before computed boundaries.

```@example
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

A committor measures the probability that a system, starting at a given point in phase space, will reach one designated region before another. Formally, for two disjoint sets A and B, the committor q(x) gives the likelihood that a trajectory initiated at x will reach B before A under the system’s dynamics.
```@example
_, Aind = find_boundary(mesh.pts, point_a, radii, density)
_, Bind = find_boundary(mesh.pts, point_b, radii, density)

q = committor(langevin_sys, mesh, Aind, Bind)

@show extrema(q)

tricontourf(Triangulation(mesh.pts', mesh.tri'), q)
```
We can also compute the reverse committor, which is the probability that a trajectory initiated at x will reach A before B under the system’s dynamics. Hence, we must reverse the drift function in the Langevin system and swap the boundaries A and B in the committor function.
```@example
function divfree1(x,y)
    f1,f2 = divfree(x,y)
    return -f1,-f2
end

langevin_sys_reverse = CriticalTransitions.Langevin(Hamiltonian, divfree1, KE, gamma, beta)

qminus = committor(langevin_sys_reverse, mesh, Bind, Aind)

@show extrema(qminus)

tricontourf(Triangulation(mesh.pts', mesh.tri'), qminus)
```



```@example

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

```@example
# probability density of reactive trajectories
mu = exp.(-beta * Hamiltonian(pts[:,1], pts[:,2])) / Z
muAB = mu .* q .* qminus

tricontourf(Triangulation(mesh.pts', mesh.tri'), muAB)
```

```@example
Rcurrent, Rrate = reactive_current(langevin_sys, mesh, q, qminus, Z)
@show Rrate
```

```@example
ARcurrent = vec(sqrt.(sum(Rcurrent.^2, dims=2)))
ARCmax = maximum(ARcurrent)

tricontourf(Triangulation(mesh.pts', mesh.tri'), ARcurrent)
```

```@example
c =ARcurrent./maxima(ARcurrent)
arrows(pts[:,1], pts[:,2], Rcurrent[:,1]./maxima(ARcurrent), Rcurrent[:,2]./maxima(ARcurrent), arrowsize = c*10, lengthscale = 0.1, arrowcolor = c, linecolor = c)
```

```@example
prob_reactive = probability_reactive(langevin_sys, mesh, q, qminus, Z)
print("Probability that a trajectory is reactive at a randomly picked time: ",prob_reactive)
```

```@example
prob_lastA = probability_last_A(langevin_sys, mesh, Amesh, qminus, Z)
print("Probability that a trajectory last visited A: ", prob_lastA)
```