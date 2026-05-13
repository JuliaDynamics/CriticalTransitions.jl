"""
$(TYPEDSIGNATURES)

Bernoulli function `B(z) = z / (exp(z) - 1)` with `B(0) = 1`. Generic in
the floating-point type of `z`; thresholds for switching between the
Taylor expansion and the direct form scale with `eps(T)` so the series
remains accurate at any precision.
"""
@inline function _bernoulli(z::T) where {T <: AbstractFloat}
    az = abs(z)
    # Below `eps^(2/3)` the linear term is already at machine precision,
    # so `B(z) ≈ 1`. Between that and `eps^(1/5)` the order-4 Taylor
    # series is accurate to ~eps. Above, use `z / expm1(z)` directly.
    near_zero = eps(T)^(T(2) / 3)
    midrange = eps(T)^(T(1) / 5)
    if az < near_zero
        return one(T)
    elseif az < midrange
        z2 = z * z
        return one(T) - z / 2 + z2 / 12 - z2 * z2 / 720
    else
        return z / expm1(z)
    end
end

"""
$(TYPEDSIGNATURES)

Impose Dirichlet rows in-place by zeroing masked rows and setting their
diagonal to one.
"""
function _enforce_dirichlet_rows!(
        A::SparseMatrixCSC{T, Int}, mask::BitVector
    ) where {T <: AbstractFloat}
    rv = rowvals(A)
    nz = nonzeros(A)
    @inbounds for col in 1:size(A, 2)
        for p in nzrange(A, col)
            mask[rv[p]] && (nz[p] = zero(T))
        end
    end
    @inbounds for n in eachindex(mask)
        mask[n] && (A[n, n] = one(T))
    end
    return A
end

"""
$(TYPEDSIGNATURES)

Solve the linear system `A q = b` subject to Dirichlet conditions
`q[mask] = values[mask]`. `source` sets the free-row right-hand side.
"""
function _solve_dirichlet(
        A::SparseMatrixCSC{T, Int},
        mask::BitVector,
        values::Vector{T},
        alg = UMFPACKFactorization();
        source::Union{Nothing, Vector{T}} = nothing,
        kwargs...,
    )::Vector{T} where {T <: AbstractFloat}
    M = copy(A)
    rhs = source === nothing ? zeros(T, size(A, 1)) : copy(source)
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

Default `alg` is `UMFPACKFactorization()` (Float64 only). Pass any
`SciMLLinearSolveAlgorithm` for an iterative solver; further kwargs
flow to `LinearSolve.solve`. For non-Float64 generators, choose an
algorithm that supports the matrix eltype (e.g. `KrylovJL_GMRES()` or
`GenericLUFactorization()`).
"""
function _invariant_density(
        G::SparseMatrixCSC{T, Int}, weights::Vector{T},
        alg = UMFPACKFactorization();
        pin_row::Int = 1, clamp_negative::Bool = true, kwargs...,
    )::Vector{T} where {T <: AbstractFloat}
    N = size(G, 1)
    length(weights) == N || throw(
        DimensionMismatch("weight vector length $(length(weights)) ≠ N = $N"),
    )
    1 <= pin_row <= N || throw(BoundsError(1:N, pin_row))

    A = sparse(transpose(G))
    rhs = zeros(T, N)
    A[pin_row, :] .= weights
    rhs[pin_row] = one(T)

    ρ = Vector{T}(solve(LinearProblem(A, rhs), alg; kwargs...).u)
    clamp_negative && clamp!(ρ, zero(T), typemax(T))
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
        G::SparseMatrixCSC{T, Int}, ρ::Vector{T}
    )::SparseMatrixCSC{T, Int} where {T <: AbstractFloat}
    N = size(G, 1)
    length(ρ) == N || throw(
        DimensionMismatch("ρ has length $(length(ρ)) but generator is $N×$N"),
    )
    # Time-reversal is only defined on a strictly positive invariant density;
    # zeros (e.g. off-basin QSD entries) need to be restricted to a support
    # mask by the caller instead of silently clamped.
    any(ρi -> ρi <= 0, ρ) && throw(
        ArgumentError(
            "ρ must be strictly positive; restrict G and ρ to the supporting " *
                "communicating class before calling `_adjoint_generator`"
        ),
    )
    Gtil = sparse(transpose(G))
    rv = rowvals(Gtil)
    nz = nonzeros(Gtil)
    diagacc = zeros(T, N)

    @inbounds for col in 1:N
        for p in nzrange(Gtil, col)
            row = rv[p]
            if row == col
                nz[p] = zero(T)
                continue
            end
            val = nz[p] * ρ[col] / ρ[row]
            nz[p] = val
            diagacc[row] -= val
        end
    end
    @inbounds for n in 1:N
        Gtil[n, n] = diagacc[n]
    end
    return Gtil
end
