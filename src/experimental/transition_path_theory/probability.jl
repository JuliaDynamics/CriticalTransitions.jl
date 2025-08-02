"""
$(TYPEDSIGNATURES)

Calculate the probability that a trajectory is reactive (transitions from A to B).

# Arguments
- `sys::LangevinSystem`: System with Hamiltonian and inverse temperature (beta)
- `mesh::Mesh`: Mesh structure containing points and triangulation
- `q`: Forward committor function values
- `qminus`: Backward committor function values
- `Z`: Partition function value

# Returns
- Probability (float) of reactive trajectories normalized by partition function
"""
function probability_reactive(sys::LangevinSystem, mesh::Mesh, q, qminus, Z)
    ham, beta = sys.Hamiltonian, sys.beta
    pts, tri = mesh.pts, mesh.tri
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
        mu = exp(-beta * ham(vmid[:, 1], vmid[:, 2])[1])
        prob += tri_area * mu * qmid * qmmid
    end
    prob /= Z
    return prob
end

function probability_reactive(sys::LangevinSystem, q::Committor)
    return probability_reactive(sys, q.mesh.mesh, q.forward, q.backward, q.Z)
end

"""
$(TYPEDSIGNATURES)

Calculate the probability that the last visited metastable state was A.

# Arguments
- `sys::LangevinSystem`: System with Hamiltonian and inverse temperature (beta)
- `mesh::Mesh`: Main mesh structure containing points and triangulation
- `Ames::Mesh`: Mesh structure for region A
- `qminus`: Backward committor function values
- `Z`: Partition function value

# Returns
- Probability (float) that the system was last in state A, normalized by partition function
"""
function probability_last_A(sys::LangevinSystem, mesh::Mesh, Ames::Mesh, qminus, Z)
    ham, beta = sys.Hamiltonian, sys.beta
    pts, tri = mesh.pts, mesh.tri
    pts_Amesh, tri_Amesh = Ames.pts, Ames.tri
    Ntri = size(tri, 1)
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
        mu = exp(-beta * ham(vmid[:, 1], vmid[:, 2])[1])
        prob += tri_area * mu * qmmid
    end
    for j in 1:Ntri_Amesh
        ind = tri_Amesh[j, :]
        verts = pts_Amesh[ind, :]
        Aux = ones(3, 3)
        Aux[2:3, :] .= transpose(verts)
        tri_area = 0.5 * abs(det(Aux))
        vmid = reshape(sum(verts; dims=1) / 3, 1, 2) # midpoint of mesh triangle
        mu = exp(-beta * ham(vmid[:, 1], vmid[:, 2])[1])
        prob += tri_area * mu
    end
    prob /= Z
    return prob
end

function probability_last_A(sys::LangevinSystem, q::Committor)
    return probability_last_A(sys, q.mesh.mesh, q.mesh.Amesh, q.backward, q.Z)
end
