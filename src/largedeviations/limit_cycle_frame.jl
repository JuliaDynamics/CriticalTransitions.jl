# LimitCycleFrame: moving frame and normal-plane coefficients along a stable limit cycle.

using DynamicalSystemsBase: TangentDynamicalSystem, current_deviations
using LinearAlgebra: I, norm, qr, schur
using StateSpaceSets: StateSpaceSet
using StaticArrays: SMatrix

"""
Return `Φ::Array{T,3}` of size `d × d × Nτ` with `Φ[:, :, k]` equal to the state
transition matrix `Φ(τ_k, 0)` of the variational equation `Φ̇ = ∂ₓb(γ(τ))·Φ` along the
periodic orbit, with `Φ(0,0) = I`.
"""
function _state_transition_matrices(sys::CoupledSDEs, points::StateSpaceSet, period::Real)
    d = length(points[1])
    Nτ = length(points)
    Δτ = period / Nτ
    Tf = eltype(points[1])
    Φs = zeros(Tf, d, d, Nτ)
    Φs[:, :, 1] = Matrix{Tf}(I, d, d)
    tands = TangentDynamicalSystem(CoupledODEs(sys); u0 = collect(points[1]))
    for k in 2:Nτ
        DynamicalSystemsBase.step!(tands, Δτ, true)
        Φs[:, :, k] = current_deviations(tands)
    end
    return Φs
end

function _qr_posR(A::AbstractMatrix)
    qf = qr(A)
    Q = Matrix(qf.Q)
    n = min(size(A)...)
    @inbounds for j in 1:n
        if qf.R[j, j] < 0
            for i in 1:size(Q, 1)
                Q[i, j] = -Q[i, j]
            end
        end
    end
    return Q
end

function _final_monodromy(sys::CoupledSDEs, points::StateSpaceSet, period::Real)
    tands = TangentDynamicalSystem(CoupledODEs(sys); u0 = collect(points[1]))
    DynamicalSystemsBase.step!(tands, period, true)
    return Matrix(current_deviations(tands))
end

"""
    LimitCycleFrame{D, T, M}

Moving frame and normal-plane coefficient matrices along a stable limit cycle Γ.

`E[:, :, k]` is expressed in the canonical (Schur) basis of the closure: column 1 is
the unit tangent, columns 2..D are arranged so that the closure of the transverse
bundle is exactly `diag(F)`. `F[j] = +1` marks a periodic transverse direction,
`F[j] = -1` marks an antiperiodic one. `S[i, j] = F[i] * F[j]` is the precomputed
sign-mask used to apply `diag(F) · X · diag(F)` as a Hadamard product.
"""
struct LimitCycleFrame{D, T, M}
    γ::StateSpaceSet{D, T}
    period::T
    E::Array{T, 3}
    Ẽ::Array{T, 3}
    M̃::Array{T, 3}
    Ã::Array{T, 3}
    F::NTuple{M, Int8}
    S::SMatrix{M, M, Int8}
end

function LimitCycleFrame(
        points::StateSpaceSet{D, T}, period::Real, sys::CoupledSDEs
    ) where {D, T}
    period_t = T(period)
    Φs = _state_transition_matrices(sys, points, period_t)
    Φ_M = _final_monodromy(sys, points, period_t)
    E, F_vec = _build_frame(sys, points, period_t, Φs, Φ_M)
    Ẽ = E[:, 2:D, :]
    M̃, Ã = _assemble_M_A(sys, points, period_t, E, Ẽ)
    M = D - 1
    F = NTuple{M, Int8}(F_vec)
    S = SMatrix{M, M, Int8}(F[i] * F[j] for i in 1:M, j in 1:M)
    return LimitCycleFrame{D, T, M}(points, period_t, E, Ẽ, M̃, Ã, F, S)
end

# Real-Schur decomposition of `R_close ∈ O(d-1)` factors as
# `R_close = Z · (U_rot_can · F_can) · Z'`, where `U_rot_can ∈ SO(d-1)` is the
# direct sum of all 2×2 rotation blocks (plus +1 entries), and `F_can` is the
# diagonal of ±1 entries. They commute (disjoint non-identity blocks), so the smear
# `exp(-(τ/T) log U_rot_can)` absorbs only the SO part. We compute the smear in
# closed form per block: no matrix log or exp is called.
function _build_frame(sys, points, period, Φs, Φ_M)
    d = length(points[1])
    Nτ = length(points)
    Tf = eltype(points[1])
    E = zeros(Tf, d, d, Nτ)
    for k in 1:Nτ
        v = drift(sys, collect(points[k]))
        E[:, 1, k] .= v ./ norm(v)
    end
    t0 = E[:, 1, 1]
    # `qr` returns Q with arbitrary column signs. Enforcing positive-diagonal R
    # makes the QR unique and the closure `Q0' * Q_T` a deterministic O(d-1)
    # change-of-basis matrix.
    Q0 = _qr_posR(hcat(t0, Matrix{Tf}(I, d, d)))[:, 2:d]
    Ẽ_raw = zeros(Tf, d, d - 1, Nτ)
    Ẽ_raw[:, :, 1] = Q0
    for k in 2:Nτ
        Yk = Φs[:, :, k] * Q0
        Yk .-= E[:, 1, k] * (E[:, 1, k]' * Yk)
        Ẽ_raw[:, :, k] = _qr_posR(Yk)[:, 1:(d - 1)]
    end
    Y_T = Φ_M * Q0
    Y_T .-= t0 * (t0' * Y_T)
    Q_T = _qr_posR(Y_T)[:, 1:(d - 1)]
    R_close = Q0' * Q_T
    sch = schur(R_close)
    Z, Tblk = sch.Z, sch.T
    F_vec = ones(Int8, d - 1)
    rot_blocks = Tuple{UnitRange{Int}, Tf}[]
    i = 1
    while i ≤ d - 1
        is_2x2 = i < d - 1 && abs(Tblk[i + 1, i]) > 1.0e-12
        if is_2x2
            θ = atan(Tblk[i + 1, i], Tblk[i, i])
            push!(rot_blocks, (i:(i + 1), θ))
            i += 2
        else
            F_vec[i] = Int8(sign(Tblk[i, i]))
            i += 1
        end
    end
    smear = Matrix{Tf}(I, d - 1, d - 1)
    for k in 1:Nτ
        α = (k - 1) / Nτ
        fill!(smear, zero(Tf))
        for j in 1:(d - 1)
            smear[j, j] = one(Tf)
        end
        for (rng, θ) in rot_blocks
            c, s = cos(α * θ), sin(α * θ)
            smear[rng[1], rng[1]] = c
            smear[rng[1], rng[2]] = s
            smear[rng[2], rng[1]] = -s
            smear[rng[2], rng[2]] = c
        end
        E[:, 2:d, k] = Ẽ_raw[:, :, k] * Z * smear
    end
    return E, F_vec
end

function _assemble_M_A(sys, points, period, E, Ẽ)
    d = length(points[1])
    Nτ = length(points)
    Δτ = period / Nτ
    a_fn = _trace_normalized_a(sys)
    drift_fn = let p = current_parameters(sys), f = dynamic_rule(sys)
        x -> f(x, p, 0.0)
    end
    M̃ = zeros(d - 1, d - 1, Nτ)
    Ã = zeros(d - 1, d - 1, Nτ)
    Erecip = zeros(d, d, Nτ)
    for k in 1:Nτ
        Erecip[:, :, k] = inv(E[:, :, k])
    end
    Ė = similar(E)
    for k in 1:Nτ
        kp = mod1(k + 1, Nτ)
        km = mod1(k - 1, Nτ)
        Ė[:, :, k] = (E[:, :, kp] - E[:, :, km]) / (2 * Δτ)
    end
    for k in 1:Nτ
        xk = collect(points[k])
        J = ForwardDiff.jacobian(drift_fn, xk)
        a_k = a_fn(xk)
        for i in 1:(d - 1), j in 1:(d - 1)
            J̃ij = dot(Erecip[i + 1, :, k], J * Ẽ[:, j, k])
            Ωij = dot(Erecip[i + 1, :, k], Ė[:, j + 1, k])
            M̃[i, j, k] = J̃ij - Ωij
            Ã[i, j, k] = dot(Erecip[i + 1, :, k], a_k * Erecip[j + 1, :, k])
        end
    end
    return M̃, Ã
end
