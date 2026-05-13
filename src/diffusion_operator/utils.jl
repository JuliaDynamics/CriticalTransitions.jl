"""
$(TYPEDSIGNATURES)

Bernoulli function `B(z) = z / (exp(z) - 1)` with `B(0) = 1`.
"""
@inline function _bernoulli(z::Float64)
    az = abs(z)
    if az < 1.0e-10
        return 1.0
    elseif az < 1.0e-2
        z2 = z * z
        return 1.0 - z / 2 + z2 / 12 - z2 * z2 / 720
    else
        return z / expm1(z)
    end
end

"""
$(TYPEDSIGNATURES)

Impose Dirichlet rows in-place by zeroing masked rows and setting their
diagonal to one.
"""
function _enforce_dirichlet_rows!(A::SparseMatrixCSC{Float64, Int}, mask::BitVector)
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
`q[mask] = values[mask]`. `source` sets the free-row right-hand side.
"""
function _solve_dirichlet(
        A::SparseMatrixCSC{Float64, Int},
        mask::BitVector,
        values::Vector{Float64},
        alg = UMFPACKFactorization();
        source::Union{Nothing, Vector{Float64}} = nothing,
        kwargs...,
    )::Vector{Float64}
    M = copy(A)
    rhs = source === nothing ? zeros(Float64, size(A, 1)) : copy(source)
    _enforce_dirichlet_rows!(M, mask)
    @inbounds for n in eachindex(mask)
        mask[n] && (rhs[n] = values[n])
    end
    return solve(LinearProblem(M, rhs), alg; kwargs...).u
end

"""
$(TYPEDSIGNATURES)

Solve for the invariant probability density `ρ` of a CTMC with generator
`G` and per-cell volume vector `weights`: `ρᵀ G = 0` augmented with the
normalisation `dot(ρ, weights) = 1`. Sign-agnostic.

Default `alg` is `UMFPACKFactorization()`. Pass any
`SciMLLinearSolveAlgorithm` for an iterative solver; further kwargs
flow to `LinearSolve.solve`.
"""
function _invariant_density(
        G::SparseMatrixCSC{Float64, Int}, weights::Vector{Float64},
        alg = UMFPACKFactorization();
        pin_row::Int = 1, clamp_negative::Bool = true, kwargs...,
    )::Vector{Float64}
    N = size(G, 1)
    length(weights) == N || throw(
        DimensionMismatch("weight vector length $(length(weights)) ≠ N = $N"),
    )
    1 <= pin_row <= N || throw(BoundsError(1:N, pin_row))

    A = sparse(transpose(G))
    rhs = zeros(Float64, N)
    A[pin_row, :] .= weights
    rhs[pin_row] = 1.0

    ρ = Vector{Float64}(solve(LinearProblem(A, rhs), alg; kwargs...).u)
    clamp_negative && clamp!(ρ, 0.0, Inf)
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
        G::SparseMatrixCSC{Float64, Int}, ρ::Vector{Float64}
    )::SparseMatrixCSC{Float64, Int}
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
