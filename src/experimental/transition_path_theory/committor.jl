
"""
$(TYPEDSIGNATURES)

Compute a local stiffness-like matrix for the Langevin dynamics on a single triangular element.

# Arguments
- `verts::AbstractMatrix{T}`: A 3×2 matrix whose rows represent the coordinates of the triangular element's vertices.

# Returns
- A 3×3 matrix of type `Float64` constituting the local contribution of this triangle to the global stiffness operator.

# Implementation Details
This function constructs an auxiliary 3×3 matrix to calculate the determinant used in local area computations. It then forms a gradient vector, reshapes it, and computes the outer product scaled by a factor involving the determinant. This local matrix is essential for capturing the diffusion-related terms in the finite element approach for the Langevin system.
"""
function stima_Langevin(verts)
    # Build a 3×3 matrix with ones and place transposed verts in rows 2 and 3
    Aux = ones(3, 3)
    Aux[2:3, :] .= transpose(verts)
    detv = abs(det(Aux))

    # Compute G (reshape to 3×1)
    G = [verts[3, 1] - verts[2, 1], verts[1, 1] - verts[3, 1], verts[2, 1] - verts[1, 1]]
    G = reshape(G, (3, 1))

    # Return local matrix
    return (0.5 / detv) * (G * transpose(G))
end

"""
$(TYPEDSIGNATURES)

Form a block of the stiffness matrix that accounts for drift or divergence-free vector fields in the Langevin PDE.

# Arguments
- `verts::AbstractMatrix{T}`: A 3×2 matrix, with each row representing a vertex of a triangle.
- `b1::T`: A scalar or array-like value representing the x-component of some drift or field at the midpoint.
- `b2::T`: A scalar or array-like value representing the y-component of some drift or field at the midpoint.

# Returns
- A 3×3 matrix that encodes how the specified vector field contributes to the finite element discretization.

# Implementation Details
The function computes adjustments for each vertex by combining the differences in y-coordinates (scaled by `b1`) with the differences in x-coordinates (scaled by `b2`). It then reshapes and replicates this row to produce a matrix whose entries incorporate movement or drift components in the element’s local contribution.
"""
function stimavbdv(verts, b1, b2)
    # Differences in y-coordinates
    bdv1 =
        b1 .*
        [verts[2, 2] - verts[3, 2], verts[3, 2] - verts[1, 2], verts[1, 2] - verts[2, 2]]

    # Differences in x-coordinates
    bdv2 =
        b2 .*
        [verts[3, 1] - verts[2, 1], verts[1, 1] - verts[3, 1], verts[2, 1] - verts[1, 1]]

    bdv_ = bdv1 .+ bdv2
    bdv = reshape(bdv_, 1, 3)

    # Replicate row 3 times
    return (1 / 6) .* vcat(bdv, bdv, bdv)
end

"""
$(TYPEDSIGNATURES)

Solve the committor equation for a two-dimensional Langevin system using finite elements.

# Arguments
- `sys::LangevinSystem`: The Langevin system containing kinetic energy, drift-free function, beta, and gamma.
- `mesh::Mesh`: The mesh containing points and triangles.
- `Aind::Vector{Int}`: Indices of mesh points corresponding to set A (Dirichlet boundary condition = 0).
- `Bind::Vector{Int}`: Indices of mesh points corresponding to set B (Dirichlet boundary condition = 1).

# Returns
- A vector of committor values of length `N`, where each entry corresponds to a mesh node.

# Implementation Details
This function assembles a global matrix `A` and right-hand-side vector `b` for the finite element discretization of the committor equation in Langevin dynamics. The code imposes Dirichlet boundary conditions on specified nodes (`Aind` set to 0, `Bind` set to 1). It computes elementwise contributions with `stima_Langevin` and `stimavbdv`, applying exponential factors involving `beta` and `gamma` to incorporate potential and damping effects. Finally, it solves the resulting linear system for the committor values on the free (non-boundary) nodes.
"""
function committor(sys::LangevinSystem, mesh::Mesh, Aind, Bind)
    KE, divfree, beta, gamma = sys.kinetic, sys.driftfree, sys.beta, sys.gamma
    pts, tri = mesh.pts, mesh.tri

    Npts = size(pts, 1)
    Ntri = size(tri, 1)
    Dir_bdry = vcat(Aind, Bind)
    free_nodes = setdiff(1:Npts, Dir_bdry)

    A = zeros(Float64, Npts, Npts)
    b = zeros(Npts, 1)
    q = zeros(Npts, 1)
    q[Bind] .= 1

    # stiffness matrix
    for j in 1:Ntri
        ind = tri[j, :]
        verts = pts[ind, :] # vertices of mesh triangle
        vmid = reshape(sum(verts; dims=1) / 3, 1, 2) # midpoint of mesh triangle
        fac = exp.(-beta * KE(vmid)[1])
        fac_gamma = gamma * fac
        fac_beta = beta * fac
        f1, f2 = divfree(vmid[1], vmid[2])

        A[ind, ind] +=
            stima_Langevin(verts) * fac_gamma - fac_beta * stimavbdv(verts, f1, f2)
    end

    # load vector
    b = b - A * q

    # solve for committor
    q[free_nodes] = A[free_nodes, free_nodes] \ b[free_nodes]
    return vec(q)
end

function _committor(sys::LangevinSystem, TPmesh::TransitionPathMesh)
    KE, divfree, beta, gamma = sys.kinetic, sys.driftfree, sys.beta, sys.gamma
    pts, tri = TPmesh.mesh.pts, TPmesh.mesh.tri

    _, Aind = find_boundary_A(TPmesh; set=:A)
    _, Bind = find_boundary_A(TPmesh; set=:B)

    return committor(sys, TPmesh.mesh, Aind, Bind)
end

struct Committor
    mesh::TransitionPathMesh
    forward::Vector{Float64}
    backward::Vector{Float64}
    Z::Float64

    function Committor(sys::LangevinSystem, mesh::TransitionPathMesh)
        forwardsys = sys
        backwardsys = remake(sys; driftfree=(x, p) -> -1 .* sys.driftfree(x, p))
        forwardmesh = mesh
        backwardmesh = reverse_AB(mesh)

        forward = _committor(forwardsys, forwardmesh)
        backward = _committor(backwardsys, backwardmesh)
        return new(mesh, forward, backward, partition_function(sys, mesh))
    end
end
function committor(sys::LangevinSystem, mesh::TransitionPathMesh)
    return Committor(sys, mesh)
end
