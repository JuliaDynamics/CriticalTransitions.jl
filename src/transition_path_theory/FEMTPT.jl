module FEMTPT

using LinearAlgebra
using SparseArrays

"""
    stima_Langevin(verts)

Compute a local stiffness-like matrix for the Langevin dynamics on a single triangular element.

# Arguments
- `verts::AbstractMatrix`: A 3×2 matrix whose rows represent the coordinates of the triangular element's vertices.

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
    stimavbdv(verts, b1, b2)

Form a block of the stiffness matrix that accounts for drift or divergence-free vector fields in the Langevin PDE.

# Arguments
- `verts::AbstractMatrix`: A 3×2 matrix, with each row representing a vertex of a triangle.
- `b1`: A scalar or array-like value representing the x-component of some drift or field at the midpoint.
- `b2`: A scalar or array-like value representing the y-component of some drift or field at the midpoint.

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
    FEM_committor_solver_Langevin(pts, tri, Aind, Bind, fpot, divfree, beta, gamma)

Solve the committor equation for a two-dimensional Langevin system using finite elements.

# Arguments
- `pts::AbstractMatrix`: An N×2 matrix containing the coordinates of the mesh points.
- `tri::AbstractMatrix`: An M×3 matrix specifying the element connectivity by indexing into `pts`.
- `Aind::Vector{Int}`: Indices of mesh points corresponding to set A (Dirichlet boundary condition = 0).
- `Bind::Vector{Int}`: Indices of mesh points corresponding to set B (Dirichlet boundary condition = 1).
- `fpot::Function`: A function that, given a point, returns the potential energy at that coordinate.
- `divfree::Function`: A function representing the divergence-free part of the drift, returning (f1, f2) at a given point.
- `beta::Real`: Inverse of temperature (related to noise intensity).
- `gamma::Real`: Damping coefficient in the Langevin dynamics.

# Returns
- A vector of committor values of length `N`, where each entry corresponds to a mesh node.

# Implementation Details
This function assembles a global matrix `A` and right-hand-side vector `b` for the finite element discretization of the committor equation in Langevin dynamics. The code imposes Dirichlet boundary conditions on specified nodes (`Aind` set to 0, `Bind` set to 1). It computes elementwise contributions with `stima_Langevin` and `stimavbdv`, applying exponential factors involving `beta` and `gamma` to incorporate potential and damping effects. Finally, it solves the resulting linear system for the committor values on the free (non-boundary) nodes.
"""

function FEM_committor_solver_Langevin(pts, tri, Aind, Bind, fpot, divfree, beta, gamma)
    # solves for committor for
    # dx = p*dt,
    # dp = (-V_x - \gamma*p)dt + \sqrt{2\beta^{-1}\gamma}dW
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
        fac = exp.(-beta * fpot(vmid)[1])
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

function invariant_pdf(pts, tri, pts_Amesh, tri_Amesh, pts_Bmesh, tri_Bmesh, fpot, beta)
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
        mu = exp(-beta * fpot(vmid)[1])
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
        mu = exp(-beta * fpot(vmid)[1])
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
        mu = exp(-beta * fpot(vmid)[1])
        Z += tri_area * mu
    end

    return Z
end

function reactive_current_transition_rate_Langevin(
    pts, tri, fpot, divfree, beta, gamma, q, qminus, Z
)
    Npts = size(pts, 1)
    Ntri = size(tri, 1)
    # find the reactive current and the transition rate
    Rcurrent = zeros(Float64, (Ntri, 2)) # reactive current at the centers of mesh triangles
    Rrate = 0.0

    for j in 1:Ntri
        ind = tri[j, :]
        verts = pts[ind, :]
        qtri = q[ind]
        qmtri = qminus[ind]

        a = [
            verts[2, 1]-verts[1, 1] verts[2, 2]-verts[1, 2]
            verts[3, 1]-verts[1, 1] verts[3, 2]-verts[1, 2]
        ]
        b = [qtri[2] - qtri[1], qtri[3] - qtri[1]]
        bm = [qmtri[2] - qmtri[1], qmtri[3] - qmtri[1]]

        g = a \ b
        gm = a \ bm

        Aux = ones(3, 3)
        Aux[2:3, :] .= transpose(verts)
        tri_area = 0.5 * abs(det(Aux))

        vmid = reshape(sum(verts; dims=1) / 3, 1, 2) # midpoint of mesh triangle
        mu = exp(-beta * fpot(vmid)[1])
        qmid = sum(qtri) / 3
        qmmid = sum(qmtri) / 3

        f1, f2 = divfree(vmid[1], vmid[2])
        Rcurrent[j, :] .= mu * qmid * qmmid * [f1, f2]
        Rcurrent[j, 1] += mu * (gamma / beta) * (qmmid * g[1] - qmid * gm[1])
        Rrate += g[1]^2 * mu * tri_area
    end

    Rrate = Rrate * gamma / (Z * beta)
    Rcurrent ./= Z

    # map reactive current on vertices
    Rcurrent_verts = zeros(Float64, (Npts, 2))
    tcount = zeros(Int, Npts) # the number of triangles adjacent to each vertex

    for j in 1:Ntri
        indt = tri[j, :]
        Rcurrent_verts[indt, :] .+= reshape(Rcurrent[j, :], 1, 2)
        tcount[indt] .+= 1
    end

    # Divide each row by the corresponding count
    for i in 1:Npts
        Rcurrent_verts[i, :] ./= tcount[i]
    end

    return Rcurrent_verts, Rrate
end

function probability_reactive_Langevin(pts, tri, fpot, beta, q, qminus, Z)
    Npts = size(pts, 1)
    Ntri = size(tri, 1)
    # find the reactive current and the transition rate
    prob = 0.0
    for j in 1:Ntri
        ind = tri[j, :]
        verts = pts[ind, :]
        qtri = q[ind]
        qmid = sum(qtri) / 3
        qmtri = qminus[ind]
        qmmid = sum(qmtri) / 3
        Aux = ones(3, 3)
        Aux[2:3, :] .= transpose(verts)
        tri_area = 0.5 * abs(det(Aux))
        vmid = reshape(sum(verts; dims=1) / 3, 1, 2) # midpoint of mesh triangle
        mu = exp(-beta * fpot(vmid)[1])
        prob += tri_area * mu * qmid * qmmid
    end
    prob /= Z
    return prob
end

function probability_last_A_Langevin(pts, tri, pts_Amesh, tri_Amesh, fpot, beta, qminus, Z)
    Npts = size(pts, 1)
    Ntri = size(tri, 1)
    Npts_Amesh = size(pts_Amesh, 1)
    Ntri_Amesh = size(tri_Amesh, 1)

    # find the reactive current and the transition rate
    prob = 0.0
    for j in 1:Ntri
        ind = tri[j, :]
        verts = pts[ind, :]
        qmtri = qminus[ind]
        qmmid = sum(qmtri) / 3
        Aux = ones(3, 3)
        Aux[2:3, :] .= transpose(verts)
        tri_area = 0.5 * abs(det(Aux))
        vmid = reshape(sum(verts; dims=1) / 3, 1, 2) # midpoint of mesh triangle
        mu = exp(-beta * fpot(vmid)[1])
        prob += tri_area * mu * qmmid
    end
    for j in 1:Ntri_Amesh
        ind = tri_Amesh[j, :]
        verts = pts_Amesh[ind, :]
        Aux = ones(3, 3)
        Aux[2:3, :] .= transpose(verts)
        tri_area = 0.5 * abs(det(Aux))
        vmid = reshape(sum(verts; dims=1) / 3, 1, 2) # midpoint of mesh triangle
        mu = exp(-beta * fpot(vmid)[1])
        prob += tri_area * mu
    end
    prob /= Z
    return prob
end

export FEM_committor_solver_Langevin,
    invariant_pdf,
    reactive_current_transition_rate_Langevin,
    probability_reactive_Langevin,
    probability_last_A_Langevin

end
