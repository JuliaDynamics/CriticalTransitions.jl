struct Mesh{T,I}
    pts::Matrix{T}
    tri::Matrix{I}
end

struct TransitionPathMesh{T,I}
    mesh::Mesh{T,I}
    Amesh::Mesh{T,I}
    Bmesh::Mesh{T,I}
    point_a::Tuple{T,T}
    point_b::Tuple{T,T}
    radii::Tuple{T,T}
    density::T
end

function find_boundary(pts, a, radii, h0)
    xc, yc = a
    rx, ry = radii
    circle = @. sqrt((pts[:, 1] - xc)^2 / rx^2 + (pts[:, 2] - yc)^2 / ry^2)
    ind = findall(circle .- 1.0 .< h0 * 1e-2)
    return length(ind), vec(ind)
end

function find_boundary_A(TPmesh::TransitionPathMesh; set=:A)
    return find_boundary(
        TPmesh.mesh.pts,
        set == :A ? TPmesh.point_a : TPmesh.point_b,
        TPmesh.radii,
        TPmesh.density,
    )
end

function reverse_AB(TPmesh::TransitionPathMesh)
    return TransitionPathMesh(
        TPmesh.mesh,
        TPmesh.Bmesh,
        TPmesh.Amesh,
        TPmesh.point_b,
        TPmesh.point_a,
        TPmesh.radii,
        TPmesh.density,
    )
end
