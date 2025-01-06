"""
$(TYPEDSIGNATURES)

Compute the invariant probability density function (PDF) for a Langevin system over a given mesh.

# Arguments
- `sys::Langevin`: The Langevin system containing the Hamiltonian and beta parameters.
- `mesh::Mesh`: The main mesh containing points and triangles.
- `Amesh::Mesh`: The mesh corresponding to region A.
- `Bmesh::Mesh`: The mesh corresponding to region B.

# Returns
- A scalar value representing the normalization constant `Z` of the invariant PDF.

# Implementation Details
This function calculates the invariant PDF by integrating the exponential of the negative Hamiltonian over the areas of the triangles in the provided meshes. The normalization constant `Z` is computed by summing the contributions from the main mesh, region A mesh, and region B mesh.
"""
function invariant_pdf(sys::Langevin, mesh::Mesh, Amesh::Mesh, Bmesh::Mesh)
    pts, tri = mesh.pts, mesh.tri
    pts_Amesh, tri_Amesh = Amesh.pts, Amesh.tri
    pts_Bmesh, tri_Bmesh = Bmesh.pts, Bmesh.tri

    beta, ham = sys.beta, sys.Hamiltonian

    Npts = size(pts, 1)
    Ntri = size(tri, 1)
    Npts_Amesh = size(pts_Amesh, 1)
    Ntri_Amesh = size(tri_Amesh, 1)
    Npts_Bmesh = size(pts_Bmesh, 1)
    Ntri_Bmesh = size(tri_Bmesh, 1)

    # find the reactive current and the transition rate
    Z = 0.0

    # Process main mesh
    for j in 1:Ntri
        ind = tri[j, :]
        verts = pts[ind, :]
        Aux = ones(3, 3)
        Aux[2:3, :] .= transpose(verts)
        tri_area = 0.5 * abs(det(Aux))
        vmid = reshape(sum(verts; dims=1) / 3, 1, 2)  # midpoint of mesh triangle
        mu = exp(-beta * ham(vmid[:, 1], vmid[:, 2])[1])
        Z += tri_area * mu
    end

    # Process A mesh
    for j in 1:Ntri_Amesh
        ind = tri_Amesh[j, :]
        verts = pts_Amesh[ind, :]
        Aux = ones(3, 3)
        Aux[2:3, :] .= transpose(verts)
        tri_area = 0.5 * abs(det(Aux))
        vmid = reshape(sum(verts; dims=1) / 3, 1, 2)
        mu = exp(-beta * ham(vmid[:, 1], vmid[:, 2])[1])
        Z += tri_area * mu
    end

    # Process B mesh
    for j in 1:Ntri_Bmesh
        ind = tri_Bmesh[j, :]
        verts = pts_Bmesh[ind, :]
        Aux = ones(3, 3)
        Aux[2:3, :] .= transpose(verts)
        tri_area = 0.5 * abs(det(Aux))
        vmid = reshape(sum(verts; dims=1) / 3, 1, 2)
        mu = exp(-beta * ham(vmid[:, 1], vmid[:, 2])[1])
        Z += tri_area * mu
    end

    return Z
end
