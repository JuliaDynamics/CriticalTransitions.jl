using Test

using CriticalTransitions
using DelaunayTriangulation
using Contour

beta = 20.0
gamma = 0.5

function Hamiltonian(x, y)
    0.5 .* y.^2 .+ 0.25 .* x.^4 .- 0.5 .* x.^2
end

function KE(x)
    0.5 .* (x[:,2].^2)
end

function divfree(x, y)
    f1 = y
    f2 = .- x.^3 .+ x
    return f1, f2
end

sys = CriticalTransitions.Langevin(Hamiltonian, divfree, KE, gamma, beta)


# pts = readdlm("test/Duffing_pts.csv", ',')
# tri = Int.(readdlm("test/Duffing_tri.csv", ','))
# mesh = CriticalTransitions.Mesh(pts, tri)

# point_a = (-1.0, 0.0)
# point_b = (1.0, 0.0)
# radii = (0.3, 0.4)
# density = 0.04

# function CriticalTransitions.find_boundary(pts, a, radii, h0)
#     xc, yc = a
#     rx, ry = radii
#     circle = @. sqrt((pts[:, 1] - xc)^2 / rx^2 + (pts[:, 2] - yc)^2 / ry^2)
#     ind = findall(circle .- 1.0 .< h0 * 1e-2)
#     return length(ind), vec(ind)
# end


# _, Aind = find_boundary(pts, point_a, radii, density)
# _, Bind = find_boundary(pts, point_b, radii, density)

function Hamiltonian(x, y)
    return 0.5 .* y .^ 2 .+ 0.25 .* x .^ 4 .- 0.5 .* x .^ 2
end

point_a = (-1.0, 0.0)
point_b = (1.0, 0.0)
radii = (0.3, 0.4)
density = 0.04

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

pts_outer = reparametrization(p_outer, density);
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

mesh = distmesh2D(dfunc, huniform, density, box, pfix)

_, Aind = find_boundary(mesh.pts, point_a, radii, density)
_, Bind = find_boundary(mesh.pts, point_b, radii, density)

q = committor(sys, mesh, Aind, Bind)

@test size(q, 1) == size(mesh.pts, 1)
@test extrema(q) == (0, 1)
