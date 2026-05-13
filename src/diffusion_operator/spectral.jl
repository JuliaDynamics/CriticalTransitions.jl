"""
$(TYPEDEF)

Dense eigensolver backend: `LinearAlgebra.eigen` on `Matrix(Q)`.

Returns every eigenpair; fine up to a few thousand cells, becomes
infeasible above that because of the dense `O(N²)` storage. Always
applicable.
"""
struct DenseEigen end

"""
$(TYPEDEF)

`KrylovKit.eigsolve` with shift-invert through LinearSolve.jl.

The workhorse for large sparse problems and for any analysis where the
target eigenvalue is not known a priori. Tuning kwargs (`tol`,
`krylovdim`, `maxiter`, …) flow to `KrylovKit.eigsolve`; pass
`inner_alg` to override the LinearSolve algorithm used for the inner
shift-invert solve.
"""
struct KrylovKitSolver end

# ---------------------------------------------------------------------
# ShiftInvertMap — apply `(A - σI)⁻¹` via a cached LinearSolve problem.
# Implements both `mul!` and the callable interface so the same type
# can drive a `mul!`-based Krylov loop and a function-as-operator
# eigensolver (KrylovKit) without duplication.
# ---------------------------------------------------------------------

struct ShiftInvertMap{T, Cache}
    n::Int
    cache::Cache
    buf::Vector{T}
end

Base.eltype(::ShiftInvertMap{T}) where {T} = T
Base.size(M::ShiftInvertMap) = (M.n, M.n)
Base.size(M::ShiftInvertMap, ::Integer) = M.n

function _apply_shift_invert!(out::AbstractVector, M::ShiftInvertMap, x::AbstractVector)
    @inbounds @simd for i in 1:M.n
        M.buf[i] = x[i]
    end
    M.cache.b = M.buf
    sol = solve!(M.cache)
    copyto!(out, sol.u)
    return out
end

function LinearAlgebra.mul!(y::AbstractVector, M::ShiftInvertMap, x::AbstractVector)
    return _apply_shift_invert!(y, M, x)
end

function (M::ShiftInvertMap)(x::AbstractVector)
    y = similar(x, M.n)
    return _apply_shift_invert!(y, M, x)
end

function _build_shift_invert(
        A::AbstractMatrix, σ::Number,
        inner_alg = UMFPACKFactorization()
    )
    n = size(A, 1)
    T = promote_type(eltype(A), typeof(σ), Float64)
    Aσ = sparse(A - σ * I)
    rhs = zeros(T, n)
    prob = LinearProblem{true}(Aσ, rhs)
    cache = init(prob, inner_alg; alias_A = false, alias_b = false)
    buf = zeros(T, n)
    return ShiftInvertMap{T, typeof(cache)}(n, cache, buf)
end

# ---------------------------------------------------------------------
# _principal_eigenpair: internal entry point used by
# `quasi_stationary_distribution`, `eigenmodes`, and the dense path of
# `stationary_distribution`. Returns `(λ, v)` with eigenvalue closest
# to `σ`.
# ---------------------------------------------------------------------

function _principal_eigenpair(A::AbstractMatrix, ::DenseEigen; σ::Number = 0.0)
    F = eigen(Matrix(A))
    i = argmin(abs.(F.values .- σ))
    return F.values[i], real.(F.vectors[:, i])
end

function _principal_eigenpair(
        A::AbstractMatrix, ::KrylovKitSolver;
        σ::Number = 0.0, v0::Union{Nothing, AbstractVector} = nothing,
        inner_alg = UMFPACKFactorization(),
        kwargs...,
    )
    n = size(A, 1)
    T = promote_type(eltype(A), typeof(σ), Float64)
    σ_eff = iszero(σ) ? T(-1.0e-12) : σ
    v_init = v0 === nothing ? normalize!(rand(T, n)) : convert(Vector{T}, v0)
    Amap = _build_shift_invert(A, σ_eff, inner_alg)

    # Eigenvalues of A closest to σ correspond to largest |μ| of the
    # shift-invert operator, where μ = 1 / (λ - σ).
    vals_μ, vecs_μ, _info = KrylovKit.eigsolve(Amap, v_init, 1, :LM; kwargs...)
    λ = σ_eff + 1 / vals_μ[1]
    return λ, real.(vecs_μ[1])
end

"""
$(TYPEDSIGNATURES)

Slowest-decaying eigenmodes of the generator `gen`. Returns `(λ, V)`
where `λ::Vector{ComplexF64}` are the eigenvalues with largest real part
(least negative — slowest decay), sorted descending, and `V` collects
the corresponding right eigenvectors of `Q` columnwise
(`Q * V[:, i] = λ[i] * V[:, i]`).

The eigenvalue closest to zero is exactly `0`, with eigenvector the
constant function (since `Q * 1 = 0`); the corresponding *left*
eigenvector is the invariant density. The remaining eigenvalues lie in
the open left half-plane, and `-1 / Re(λ_k)` gives the metastable
timescale of the `k`-th mode. The slow non-trivial right eigenvectors
are the canonical reaction coordinates / metastable-state indicators
used in Markov state model decomposition (e.g. PCCA+).

# Algorithm (`alg`)

- [`KrylovKitSolver()`](@ref) (default) — shift-invert Arnoldi via
  `KrylovKit.eigsolve` near `σ = 0`. The right choice for large sparse
  generators. Kwargs flow to `KrylovKit.eigsolve`; pass `inner_alg`
  (a `SciMLLinearSolveAlgorithm`, defaults to `UMFPACKFactorization()`)
  to override the LinearSolve algorithm for the inner shift-invert
  solve.
- [`DenseEigen()`](@ref) — `LinearAlgebra.eigen` on `Matrix(gen.Q)`,
  returning *all* eigenpairs sliced to the first `k`. Fine up to
  ~5000 cells.

Passing a `LinearSolve.SciMLLinearSolveAlgorithm` is rejected: the
target eigenvalues are unknown a priori, so a single linear solve is
insufficient.
"""
function eigenmodes(
        gen::DiffusionGenerator, k::Integer = 10,
        alg::Union{DenseEigen, KrylovKitSolver} = KrylovKitSolver(); kwargs...
    )
    N = size(gen.Q, 1)
    k >= 1 || throw(ArgumentError("k must be ≥ 1"))
    k = min(k, N)
    return _eigenmodes(gen.Q, k, alg; kwargs...)
end

function eigenmodes(
        ::DiffusionGenerator, _k::Integer,
        alg; kwargs...
    )
    throw(
        ArgumentError(
            "eigenmodes: a `SciMLLinearSolveAlgorithm` is not a valid backend " *
                "because the target eigenvalues are unknown a priori. " *
                "Use `KrylovKitSolver()` (default) or `DenseEigen()`."
        )
    )
end

function _eigenmodes(Q::SparseMatrixCSC, k::Integer, ::DenseEigen; kwargs...)
    F = eigen(Matrix(Q))
    perm = sortperm(F.values; by = real, rev = true)
    inds = perm[1:k]
    return F.values[inds], F.vectors[:, inds]
end

function _eigenmodes(
        Q::SparseMatrixCSC, k::Integer, ::KrylovKitSolver;
        inner_alg = UMFPACKFactorization(),
        v0::Union{Nothing, AbstractVector} = nothing,
        kwargs...,
    )
    n = size(Q, 1)
    σ = -1.0e-12  # tiny shift away from the exact 0 eigenvalue
    v_init = v0 === nothing ? normalize!(rand(Float64, n)) : convert(Vector{Float64}, v0)
    Amap = _build_shift_invert(Q, σ, inner_alg)
    vals_μ, vecs_μ, _info = KrylovKit.eigsolve(Amap, v_init, k, :LM; kwargs...)
    # μ = 1 / (λ - σ); eigenvalues of Q closest to σ ≈ 0 correspond to
    # largest |μ|, which is what `:LM` selects.
    vals = ComplexF64[σ + 1 / μ for μ in vals_μ[1:k]]
    V = hcat(vecs_μ[1:k]...)
    perm = sortperm(vals; by = real, rev = true)
    return vals[perm], V[:, perm]
end
