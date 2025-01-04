
function probability_reactive(sys::Langevin, mesh::Mesh, q, qminus, Z)
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

function probability_last_A(sys::Langevin, mesh::Mesh, Ames::Mesh, qminus, Z)
    ham, beta = sys.Hamiltonian, sys.beta
    pts, tri = mesh.pts, mesh.tri
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
