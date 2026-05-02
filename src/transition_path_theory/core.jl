# =====================================================================
# Internal numerical helpers shared across [`DiffusionChain`](@ref) and
# [`TPT`](@ref). Operators are passed in CTMC convention: positive
# off-diagonals (transition rates), negative diagonal, row sums zero.
# Most helpers are sign-agnostic — they work on any sparse matrix that
# annihilates the constant function — so they can also be used directly
# on the M-matrix view (`-chain.generator`) when convenient. Only
# [`_reactive_rate`](@ref) depends on the sign convention; it expects
# the CTMC form.
# =====================================================================

"""
$(TYPEDSIGNATURES)

Bernoulli function `B(z) = z / (eᶻ - 1)` with `B(0) = 1`. Uses a Taylor
expansion for `|z| < 1e-2` to preserve precision near the origin.

Underpins the Scharfetter–Gummel exponential-fitting finite-volume scheme:
it interpolates between centered-difference (Pe → 0) and upwind (Pe → ∞)
face stencils.
"""
@inline function _bernoulli(z::Float64)
    az = abs(z)
    if az < 1e-10
        return 1.0
    elseif az < 1e-2
        z2 = z * z
        return 1.0 - z / 2 + z2 / 12 - z2 * z2 / 720
    else
        return z / expm1(z)
    end
end

"""
$(TYPEDSIGNATURES)

In-place: zero all structural non-zeros in rows where `mask` is `true`,
then set the diagonal of those rows to `1`. Used to impose Dirichlet
boundary conditions on a sparse linear system without disturbing the
sparsity pattern of the free rows.
"""
function _enforce_dirichlet_rows!(A::SparseMatrixCSC{Float64,Int}, mask::BitVector)
    rv = rowvals(A)
    nz = nonzeros(A)
    @inbounds for col in 1:size(A, 2)
        for p in nzrange(A, col)
            mask[rv[p]] && (nz[p] = 0.0)
        end
    end
    @inbounds for n in eachindex(mask)
        mask[n] && (A[n, n] = 1.0)
    end
    return A
end

"""
$(TYPEDSIGNATURES)

Solve the linear system `A q = b` subject to Dirichlet conditions
`q[mask] = values[mask]`, returning the full `q::Vector{Float64}` of length
`size(A, 1)`. Pass `source` to set the right-hand side on free rows;
defaults to zero (i.e. solves the homogeneous problem `A q = 0` with
Dirichlet BCs). Sign-agnostic — works for both CTMC generators and PDE
M-matrices.
"""
function _solve_dirichlet(
    A::SparseMatrixCSC{Float64,Int},
    mask::BitVector,
    values::Vector{Float64};
    source::Union{Nothing,Vector{Float64}}=nothing,
)::Vector{Float64}
    M = copy(A)
    rhs = source === nothing ? zeros(Float64, size(A, 1)) : copy(source)
    _enforce_dirichlet_rows!(M, mask)
    @inbounds for n in eachindex(mask)
        mask[n] && (rhs[n] = values[n])
    end
    return M \ rhs
end

"""
$(TYPEDSIGNATURES)

Solve for the invariant probability density `ρ` of a CTMC with generator
`G` and per-cell volume vector `weights`: `ρᵀ G = 0` augmented with the
normalisation `dot(ρ, weights) = 1`. Sign-agnostic.
"""
function _invariant_density(
    G::SparseMatrixCSC{Float64,Int}, weights::Vector{Float64}
)::Vector{Float64}
    N = size(G, 1)
    length(weights) == N || throw(
        DimensionMismatch("weight vector length $(length(weights)) ≠ N = $N"),
    )

    A = sparse(transpose(G))
    rhs = zeros(Float64, N)
    A[1, :] .= weights
    rhs[1] = 1.0

    ρ = A \ rhs
    clamp!(ρ, 0.0, Inf)
    tot = dot(ρ, weights)
    tot > 0 && (@. ρ = ρ / tot)
    return ρ
end

"""
$(TYPEDSIGNATURES)

Build the time-reversed CTMC generator from `G` and the invariant density `ρ`:

    G̃[i, j] = G[j, i] · ρ[j] / ρ[i]   for i ≠ j ,

with diagonal `G̃[i, i] = -∑_{j ≠ i} G̃[i, j]`. Output preserves the
sparsity pattern of `Gᵀ` and the sign convention of `G`.
"""
function _adjoint_generator(
    G::SparseMatrixCSC{Float64,Int}, ρ::Vector{Float64}
)::SparseMatrixCSC{Float64,Int}
    N = size(G, 1)
    Gtil = sparse(transpose(G))
    rv = rowvals(Gtil)
    nz = nonzeros(Gtil)
    diagacc = zeros(Float64, N)
    ρ_safe = max.(ρ, eps(Float64))

    @inbounds for col in 1:N
        for p in nzrange(Gtil, col)
            row = rv[p]
            if row == col
                nz[p] = 0.0
                continue
            end
            val = nz[p] * ρ_safe[col] / ρ_safe[row]
            nz[p] = val
            diagacc[row] -= val
        end
    end
    @inbounds for n in 1:N
        Gtil[n, n] = diagacc[n]
    end
    return Gtil
end

"""
$(TYPEDSIGNATURES)

A → B reactive transition rate from a CTMC generator `G` (positive
off-diagonals = transition rates):

    k_{AB} = ∑_{i ∈ A, j ∉ A} ρ[i] · v[i] · q⁻[i] · G[i, j] · q⁺[j] .

Computed in `O(nnz(G))` by walking sparse columns once.
"""
function _reactive_rate(
    G::SparseMatrixCSC{Float64,Int},
    ρ::Vector{Float64},
    qplus::Vector{Float64},
    qminus::Vector{Float64},
    weights::Vector{Float64},
    A_mask::BitVector,
)::Float64
    rv = rowvals(G)
    nz = nonzeros(G)
    rate = 0.0
    @inbounds for col in 1:size(G, 2)
        A_mask[col] && continue
        qpc = qplus[col]
        for p in nzrange(G, col)
            row = rv[p]
            row != col || continue
            A_mask[row] || continue
            rate += ρ[row] * weights[row] * qminus[row] * nz[p] * qpc
        end
    end
    return rate
end
