module DistMesh2D

using CriticalTransitions
using DelaunayTriangulation
using Interpolations
using LinearAlgebra

using CriticalTransitions: Mesh

export distmesh2D, dellipse, ddiff
export get_ellipse, reparametrization
export find_boundary, huniform, dunion

"""
    huniform(p::AbstractMatrix)

Return a uniform length scale for each point in p, suitable for a constant mesh density.

Arguments:
- p: An N×2 matrix of point coordinates.

Returns:
- Vector of uniform ones.

"""
function CriticalTransitions.huniform(p)
    m, _ = size(p)
    return ones(m, 1)
end

function CriticalTransitions.ddiff(d1, d2)
    return maximum([d1 -d2]; dims=2)
end

function CriticalTransitions.dellipse(p, a, radii)
    xc, yc = a
    rx, ry = radii
    return sqrt.(((p[:, 1] .- xc) .^ 2) / rx^2 + ((p[:, 2] .- yc) .^ 2) / ry^2) .- 1
end

function CriticalTransitions.dunion(d1, d2)
    return [minimum([d1[i], d2[i]]) for i in 1:length(d1)]
    # minimum.([d1 d2], dims=2)
end

function triarea(pts, tri)
    d12 = pts[tri[:, 2], :] - pts[tri[:, 1], :]
    d13 = pts[tri[:, 3], :] - pts[tri[:, 1], :]
    return d12[:, 1] .* d13[:, 2] - d12[:, 2] .* d13[:, 1]
end

function find_boundary(pts, a, radii, h0)
    xc, yc = a
    rx, ry = radii
    circle = @. sqrt((pts[:, 1] - xc)^2 / rx^2 + (pts[:, 2] - yc)^2 / ry^2)
    ind = findall(circle .- 1.0 .< h0 * 1e-2)
    return length(ind), vec(ind)
end

"""
    fixmesh(pts::AbstractMatrix, tri::AbstractMatrix)

Remove duplicate nodes, reorder triangle connectivity for negative areas,
and discard triangles with very small areas. Also update node indices.

Arguments:
- pts: N×2 matrix of node coordinates
- tri: M×3 matrix of triangle indices

Returns:
- Updated pts, tri without duplicate/unused nodes
"""
function fixmesh(pts, tri)
    TOL = 1.0e-10
    # Remove repeated nodes
    unique_pts, idx = unique_with_indices(pts)
    tri = reshape(idx[tri], size(tri))

    # Compute areas and reorder negative area triangles
    A = triarea(pts, tri)
    idx_tri_reorder = findall(A .< 0)
    if !isempty(idx_tri_reorder)
        tmp = tri[idx_tri_reorder, 1]
        tri[idx_tri_reorder, 1] = tri[idx_tri_reorder, 2]
        tri[idx_tri_reorder, 2] = tmp
    end

    # Remove small area triangles
    idx_keep = findall(abs.(A) .> TOL * norm(A, Inf))
    tri = tri[idx_keep, :]

    # Remove unused nodes
    t_col = reshape(tri, :)
    used_nodes = unique(t_col)
    old_to_new = zeros(Int, maximum(used_nodes))
    old_to_new[used_nodes] = 1:length(used_nodes)
    pts = pts[used_nodes, :]
    tri = reshape(old_to_new[t_col], size(tri))

    return pts, tri
end

"""
    distmesh2D(fd::Function, fh::Function, h0::Real, bbox::AbstractVector, pfix::AbstractMatrix)

Generate a 2D triangular mesh in a region defined by the distance function `fd`.
The scaled edge length function `fh` can control mesh density. The bounding box
`bbox` is used for initial point placement, and any fixed points `pfix` are added
before mesh generation.

Arguments:
- fd: Distance function d(x, y) < 0 inside region, = 0 on boundary
- fh: Scaled edge length function h(x, y)
- h0: Initial edge length
- bbox: Vector [xmin, xmax, ymin, ymax]
- pfix: Fixed node positions (each row is a point)

Returns:
- pts: Matrix of node positions
- tri: Triangle indices (connectivity)
"""
function CriticalTransitions.distmesh2D(fd, fh, h0, bbox, pfix)
    # Parameters
    dptol = 0.001
    ttol = 0.1
    Fscale = 1.2
    deltat = 0.2
    geps = 0.001 * h0
    deps = sqrt(eps(Float64)) * h0
    MAXcount = 5000
    densityctrlfreq = 30
    jshow = 200

    # Create initial distribution of points
    ax = range(bbox[1], bbox[2]; step=h0)
    ay = range(bbox[3], bbox[4]; step=h0 * sqrt(3) * 0.5)
    x = [i for i in ax, j in ay]
    y = [j for i in ax, j in ay]
    nx, ny = size(x)

    # Shift odd rows
    for i in 2:2:nx
        x[i, :] .+= h0 * 0.5
    end

    pts = hcat(vec(x), vec(y))

    # Remove points outside region
    mask = vec(fd(pts) .> geps)
    pts = pts[.!mask, :]

    # Add fixed points
    if !isempty(pfix)
        pfix = unique(pfix; dims=1)
        nfix = size(pfix, 1)
        pts = vcat(pfix, pts)
    else
        nfix = 0
    end

    count = 0
    pts_old = fill(Inf, size(pts))

    while count < MAXcount
        count += 1

        # Retriangulation if points have moved too much
        if maximum(sqrt.(sum((pts .- pts_old) .^ 2; dims=2))) / h0 > ttol
            pts_old = copy(pts)

            # Delaunay triangulation
            tri = DelaunayTriangulation.triangulate(pts')
            tri = Matrix(
                reduce(hcat, [[triangle_vertices(T)...] for T in each_solid_triangle(tri)])'
            )

            # Get centroids and remove triangles outside
            pmid = (pts[tri[:, 1], :] .+ pts[tri[:, 2], :] .+ pts[tri[:, 3], :]) / 3
            tri = tri[vec(fd(pmid) .< -geps), :]

            # Create bars from triangles
            bars = vcat(tri[:, 1:2], tri[:, 2:3], hcat(tri[:, 3], tri[:, 1]))
            bars = unique(sort(bars; dims=2); dims=1)

            # Calculate bar vectors and lengths
            barvec = pts[bars[:, 1], :] .- pts[bars[:, 2], :]
            L = sqrt.(sum(barvec .^ 2; dims=2))
            L0 = fh((pts[bars[:, 1], :] .+ pts[bars[:, 2], :]) / 2)
            L0 *= Fscale * sqrt(sum(L .^ 2) / sum(L0 .^ 2))

            # Density control
            if count % densityctrlfreq == 0
                L0_too_large = findall(vec(L0) .> 2vec(L))
                if !isempty(L0_too_large)
                    qremove = unique(reshape(bars[L0_too_large, :], :))
                    qremove = setdiff(qremove, 1:nfix)
                    pts = pts[setdiff(1:size(pts, 1), qremove), :]
                    pts_old = fill(Inf, size(pts))
                    continue
                end
            end

            # Force vector calculation
            F = @. max(L0 - L, 0)
            Fvec = F ./ L .* barvec

            # Sum up forces at nodes
            Ftot = zeros(size(pts))
            for i in 1:size(bars, 1)
                Ftot[bars[i, 1], :] .+= Fvec[i, :]
                Ftot[bars[i, 2], :] .-= Fvec[i, :]
            end

            # Zero force for fixed points
            Ftot[1:nfix, :] .= 0

            # Update points
            pts .+= deltat * Ftot

            # Project points outside to boundary
            d = fd(pts)
            ix = getindex.(findall(d .> 0), 1)

            if !isempty(ix)
                dgradx = (fd(pts[ix, :] .+ [deps 0]) .- d[ix]) / deps
                dgrady = (fd(pts[ix, :] .+ [0 deps]) .- d[ix]) / deps
                pts[ix, :] -= d[ix] .* [dgradx dgrady] ./ (dgradx .^ 2 .+ dgrady .^ 2)
            end

            # Check convergence
            d = fd(pts)
            ix = findall(d .< -geps)
            displacement = maximum(sqrt.(sum((deltat * Ftot[ix, :]) .^ 2; dims=2))) / h0

            if displacement < dptol
                break
            end

            # Show progress
            if count % jshow == 0
                println("Iteration $count: displacement = $displacement")
            end
        end
    end

    # Final cleanup
    pts, tri = fixmesh(pts, tri)

    return Mesh(pts, tri)
end

# Helper functions
function unique_with_indices(A)
    uniqueA = unique(A; dims=1)
    indices = [findfirst(r -> all(r .== row), eachrow(uniqueA)) for row in eachrow(A)]
    return uniqueA, indices
end

function CriticalTransitions.get_ellipse(a, raddii, n)
    t = range(0, 2π; length=n + 1)
    pts = zeros(n, 2)
    xc, yc = a
    rx, ry = raddii
    for i in 1:n
        pts[i, 1] = xc + rx * cos(t[i])
        pts[i, 2] = yc + ry * sin(t[i])
    end
    return pts
end

function CriticalTransitions.reparametrization(path, h)
    dp = path .- circshift(path, (1, 0))
    dp[1, :] .= 0
    dl = sqrt.(sum(dp .^ 2; dims=2))
    lp = vec(cumsum(dl; dims=1))
    total_len = lp[end]
    lp ./= total_len
    npath = round(Int, total_len / h)
    g1 = range(0, 1; length=npath)
    itp_x = LinearInterpolation(lp, path[:, 1])
    itp_y = LinearInterpolation(lp, path[:, 2])
    path_x = itp_x.(g1)
    path_y = itp_y.(g1)
    return [path_x path_y]
end
end # module
