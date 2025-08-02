using Test

using CriticalTransitions
using DelaunayTriangulation
using Contour

using CriticalTransitions:
    dellipse, get_ellipse, ddiff, dunion, reparameterization, distmesh2D

function Hamiltonian(x, y)
    return 0.5 .* y .^ 2 .+ 0.25 .* x .^ 4 .- 0.5 .* x .^ 2
end

point_a = (-1.0, 0.0)
point_b = (1.0, 0.0)
radii = (0.3, 0.4)
density = 0.1

Na = round(Int, π * sum(radii) / density) # the number of points on the A-circle
Nb = Na

ptsA = get_ellipse(point_a, radii, Na)
ptsB = get_ellipse(point_b, radii, Na);

Hbdry = 0.5
nx, ny = 41, 41
nxy = nx * ny
xmin, xmax = -2.0, 2.0
ymin, ymax = -2.0, 2.0

x1 = range(xmin, xmax; length=nx)
y1 = range(ymin, ymax; length=ny)
x_grid = [xx for yy in y1, xx in x1]
y_grid = [yy for yy in y1, xx in x1]

Hgrid = Hamiltonian(x_grid, y_grid)

cont = Contour.contour(x1, y1, Hgrid, Hbdry)
yc, xc = coordinates(lines(cont)[1])
p_outer = [xc yc]

pts_outer = reparameterization(p_outer, density);
Nouter = size(pts_outer, 1)
Nfix = Na + Nb + Nouter

box = [xmin, xmax, ymin, ymax]
pfix = zeros(Nfix, 2)
pfix[1:Na, :] .= ptsA
pfix[(Na + 1):(Na + Nb), :] .= ptsB
pfix[(Na + Nb + 1):Nfix, :] .= pts_outer

function dfunc(p)
    d0 = Hamiltonian(p[:, 1], p[:, 2])
    dA = dellipse(p, point_a, radii)
    dB = dellipse(p, point_b, radii)
    d = ddiff(d0 .- Hbdry, dunion(dA, dB))
    return d
end

mesh = distmesh2D(dfunc, CriticalTransitions.huniform, density, box, pfix)

@test size(mesh.pts, 1) == 872
@test size(mesh.tri, 1) == 1586 || size(mesh.tri, 1) == 1587

function triarea(mesh)
    pts, tri = mesh.pts, mesh.tri
    d12 = pts[tri[:, 2], :] - pts[tri[:, 1], :]
    d13 = pts[tri[:, 3], :] - pts[tri[:, 1], :]
    return d12[:, 1] .* d13[:, 2] - d12[:, 2] .* d13[:, 1]
end

@test all(triarea(mesh) .< density)
