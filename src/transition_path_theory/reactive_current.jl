"""
$(TYPEDSIGNATURES)

Calculate the reactive current and transition rate for a Langevin system.

The reactive current represents the probability flux of reactive trajectories
that successfully transition from state A to state B.

# Arguments
- `sys::Langevin`: System containing Hamiltonian, drift-free component, inverse temperature (beta), and friction coefficient (gamma)
- `mesh::Mesh`: Mesh structure containing points and triangulation
- `q`: Forward committor function values
- `qminus`: Backward committor function values
- `Z`: Partition function value

# Returns
A tuple containing:
- `Rcurrent_verts`: Array of reactive current vectors at mesh vertices
- `Rrate`: Transition rate between states A and B
"""
function reactive_current(sys::Langevin, mesh::Mesh, q, qminus, Z)
    ham, divfree, beta, gamma = sys.Hamiltonian, sys.driftfree, sys.beta, sys.gamma
    pts, tri = mesh.pts, mesh.tri
    Npts = size(pts, 1)
    Ntri = size(tri, 1)
    # find the reactive current and the transition rate
    # reactive current at the centers of mesh triangles
    Rcurrent = zeros(Float64, (Ntri, 2))
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
        mu = exp(-beta * ham(vmid[:, 1], vmid[:, 2])[1])
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

function reactive_current(sys::Langevin, q::Committor)
    return reactive_current(sys::Langevin, q.mesh.mesh::Mesh, q.forward, q.backward, q.Z)
end

function reactive_rate(sys::Langevin, mesh::Mesh, q, Z)
    ham, divfree, beta, gamma = sys.Hamiltonian, sys.driftfree, sys.beta, sys.gamma
    pts, tri = mesh.pts, mesh.tri
    Npts = size(pts, 1)
    Ntri = size(tri, 1)
    # find the reactive current and the transition rate
    # reactive current at the centers of mesh triangles
    Rcurrent = zeros(Float64, (Ntri, 2))
    Rrate = 0.0

    for j in 1:Ntri
        ind = tri[j, :]
        verts = pts[ind, :]
        qtri = q[ind]

        a = [
            verts[2, 1]-verts[1, 1] verts[2, 2]-verts[1, 2]
            verts[3, 1]-verts[1, 1] verts[3, 2]-verts[1, 2]
        ]
        b = [qtri[2] - qtri[1], qtri[3] - qtri[1]]

        g = a \ b

        Aux = ones(3, 3)
        Aux[2:3, :] .= transpose(verts)
        tri_area = 0.5 * abs(det(Aux))

        mu = exp(-beta * ham(vmid[:, 1], vmid[:, 2])[1])
        Rrate += g[1]^2 * mu * tri_area
    end

    Rrate = Rrate * gamma / (Z * beta)

    return Rrate
end

function reactive_rate(sys::Langevin, q::Committor)
    return reactive_rate(sys::Langevin, q.mesh.mesh::Mesh, q.forward, q.Z)
end
