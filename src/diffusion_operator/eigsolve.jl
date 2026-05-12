# =====================================================================
# Algorithms for the invariant density `Qᵀ ρ = 0` and, more generally,
# the principal eigenpair of a sparse generator-like matrix.
#
# Type hierarchy:
#
#     StationaryAlgorithm                       (abstract)
#     ├── LUPinned                              direct LU pin
#     └── EigensolveAlgorithm                   (abstract; also solves a general
#         ├── ArnoldiSI                              principal-eigenpair problem)
#         ├── ShiftInvertPower
#         └── KrylovKitArnoldi                  (defined here; method in extension)
#
# Every concrete struct accepts the same two LinearSolve.jl passthroughs
# via keyword constructor:
#
#     linsolve         — a SciMLLinearSolveAlgorithm or `nothing`
#                        (defaults to UMFPACK on sparse).
#     linsolve_kwargs  — `NamedTuple` of options forwarded to `init`/`solve`,
#                        e.g. `(; abstol = 1e-14, reltol = 1e-14)`.
#
# Functions consuming these:
#
#     stationary_distribution(gen)                                 → LUPinned()
#     stationary_distribution(gen, alg::StationaryAlgorithm)       → dispatch
#     quasi_stationary_distribution(gen, basin)                    → ArnoldiSI()
#     quasi_stationary_distribution(gen, basin, alg::EigensolveAlgorithm)
#     principal_eigenpair(A, alg::EigensolveAlgorithm; σ, v0)
#
# =====================================================================

"""
$(TYPEDEF)

Abstract supertype of algorithms that solve `Qᵀ ρ = 0` for the invariant
density of a discrete generator (the principal eigenvector at `λ = 0`).
Concrete subtype: [`LUPinned`](@ref). All [`EigensolveAlgorithm`](@ref)
subtypes are also `StationaryAlgorithm`s because the same problem can be
solved as a shifted principal-eigenpair problem at `σ = 0`.

[`stationary_distribution`](@ref) accepts any `StationaryAlgorithm`.
[`quasi_stationary_distribution`](@ref) requires an
[`EigensolveAlgorithm`](@ref) because the basin sub-generator has a
non-zero exit rate.
"""
abstract type StationaryAlgorithm end

"""
$(TYPEDEF)

Abstract supertype of algorithms that compute a principal eigenpair
(the eigenvalue closest to a shift `σ` and the corresponding right
eigenvector) of a sparse matrix. Concrete subtypes: [`ArnoldiSI`](@ref),
[`ShiftInvertPower`](@ref), `KrylovKitArnoldi` (available when
`KrylovKit` is loaded).

`EigensolveAlgorithm <: StationaryAlgorithm`: every eigensolver
implicitly solves the invariant-density problem at `σ = 0`.
"""
abstract type EigensolveAlgorithm <: StationaryAlgorithm end

# ---------------------------------------------------------------------
# LUPinned — direct LU pin
# ---------------------------------------------------------------------

"""
$(TYPEDEF)

Direct LU pin for the invariant density: replaces row `pin_row` of `Qᵀ`
with the normalisation constraint `∫ ρ dV = 1` and solves the augmented
sparse linear system with LinearSolve.jl.

Fast and exact for non-degenerate generators (a single near-zero
eigenvalue of `Qᵀ`). At very small noise the generator can have
multiple near-zero eigenvalues within machine precision (multi-basin
metastability), in which case the answer becomes solver- and
pin-dependent. The post-solve diagnostic in
[`stationary_distribution`](@ref) detects this and warns.

# Fields
$(TYPEDFIELDS)

# Example
```julia
# Default: sparse direct LU
stationary_distribution(gen, LUPinned())

# Iterative inner solve with tight tolerance
using LinearSolve
alg = LUPinned(; linsolve = KrylovJL_GMRES(),
                 linsolve_kwargs = (; abstol = 1e-14, reltol = 1e-14))
stationary_distribution(gen, alg)
```
"""
struct LUPinned{A, K <: NamedTuple} <: StationaryAlgorithm
    "LinearSolve.jl algorithm for the inner solve. `nothing` uses sparse direct LU."
    linsolve::A
    "Options forwarded to `LinearSolve.init` / `solve!` (e.g. `(; abstol = 1e-14)`)."
    linsolve_kwargs::K
    "Row of `Qᵀ` replaced by the normalisation constraint. Default `1`."
    pin_row::Int
end

LUPinned(; linsolve = nothing, linsolve_kwargs::NamedTuple = (;), pin_row::Int = 1) =
    LUPinned(linsolve, linsolve_kwargs, pin_row)

# ---------------------------------------------------------------------
# ShiftInvertMap — a `mul!`-compatible wrapper that applies `(A - σI)⁻¹`
# via LinearSolve.jl. Driven by Arnoldi (`ArnoldiSI`) or by the simple
# inverse-power loop (`ShiftInvertPower`).
# ---------------------------------------------------------------------

struct ShiftInvertMap{T, Cache}
    n::Int
    cache::Cache
    buf::Vector{T}
end

Base.eltype(::ShiftInvertMap{T}) where {T} = T
Base.size(M::ShiftInvertMap) = (M.n, M.n)
Base.size(M::ShiftInvertMap, ::Integer) = M.n

function LinearAlgebra.mul!(y::AbstractVector, M::ShiftInvertMap, x::AbstractVector)
    @inbounds @simd for i in 1:M.n
        M.buf[i] = x[i]
    end
    M.cache.b = M.buf
    sol = solve!(M.cache)
    copyto!(y, sol.u)
    return y
end

function _build_shift_invert(A::AbstractMatrix, σ::Number, linsolve, linsolve_kwargs::NamedTuple)
    n = size(A, 1)
    T = promote_type(eltype(A), typeof(σ), Float64)
    Aσ = sparse(A - σ * I)
    rhs = zeros(T, n)
    prob = LinearProblem{true}(Aσ, rhs)
    cache = linsolve === nothing ?
        init(prob; alias_A = false, alias_b = false, linsolve_kwargs...) :
        init(prob, linsolve; alias_A = false, alias_b = false, linsolve_kwargs...)
    buf = zeros(T, n)
    return ShiftInvertMap{T, typeof(cache)}(n, cache, buf)
end

# ---------------------------------------------------------------------
# Custom Arnoldi (no eigensolver dep), adapted from QuantumToolbox.jl.
# Builds an `m`-step Arnoldi factorisation
#
#       A V_m = V_m H_m + β v_{m+1} eₘᵀ ,
#
# does an in-place Schur reorder of `H_m` to bring the wanted Ritz
# values to the leading block, and restarts. Returns when the leading
# `k` residuals are below `tol` or after `maxiter` restarts.
# ---------------------------------------------------------------------

function _arnoldi_init!(A, b::AbstractVector{T},
        V::AbstractMatrix{T}, H::AbstractMatrix{T}) where {T <: Number}
    v₁ = view(V, :, 1)
    v₂ = view(V, :, 2)
    v₁ .= b
    normalize!(v₁)

    mul!(v₂, A, v₁)
    H[1, 1] = dot(v₁, v₂)
    axpy!(-H[1, 1], v₁, v₂)
    H[2, 1] = norm(v₂)
    v₂ ./= H[2, 1]
    return nothing
end

function _arnoldi_step!(A, V::AbstractMatrix{T}, H::AbstractMatrix{T},
        i::Integer) where {T <: Number}
    vᵢ  = view(V, :, i)
    vᵢ₊ = view(V, :, i + 1)
    mul!(vᵢ₊, A, vᵢ)
    for j in 1:i
        vⱼ = view(V, :, j)
        H[j, i] = dot(vⱼ, vᵢ₊)
        axpy!(-H[j, i], vⱼ, vᵢ₊)
    end
    β = norm(vᵢ₊)
    H[i + 1, i] = β
    if β > 0
        vᵢ₊ ./= β
    end
    return β
end

function _update_schur_eigs!(Hₘ, Uₘ, Uₘᵥ, f, k, β, sorted_vals, sortby, rev)
    copyto!(Uₘ, Hₘ)
    F = schur!(Uₘ)

    values = F.values
    sortperm!(sorted_vals, values, by = sortby, rev = rev)

    select = fill(false, length(values))
    @inbounds for j in 1:k
        select[sorted_vals[j]] = true
    end
    ordschur!(F, select)

    copyto!(Hₘ, F.T)
    copyto!(Uₘ, F.Z)
    mul!(f, Uₘᵥ, β)
    return nothing
end

# Pure-Julia right eigenvectors from upper-triangular Schur form.
function _schur_right_eigenvectors(Tₘ::AbstractMatrix, k::Integer)
    n = size(Tₘ, 1)
    vecs = zeros(eltype(Tₘ), n, k)

    @inbounds for col in 1:k
        vec = view(vecs, :, col)
        vec[col] = one(eltype(Tₘ))
        λ = Tₘ[col, col]
        for row in (col - 1):-1:1
            acc = zero(eltype(Tₘ))
            for inner in (row + 1):col
                acc += Tₘ[row, inner] * vec[inner]
            end
            denom = Tₘ[row, row] - λ
            tol = eps(typeof(abs(λ))) * max(one(typeof(abs(λ))), abs(λ))
            vec[row] = if abs(denom) <= tol
                zero(eltype(Tₘ))
            else
                -acc / denom
            end
        end
        normalize!(vec)
    end
    return vecs
end

function _eigsolve_arnoldi(A, v₀::AbstractVector{T}, k::Int, m::Int;
        tol::Real = 1.0e-10, maxiter::Int = 200,
        sortby::Function = abs, rev::Bool = true) where {T <: Number}

    n = size(A, 2)
    V = similar(v₀, n, m + 1)
    H = zeros(T, m + 1, m)

    _arnoldi_init!(A, v₀, V, H)
    @inbounds for i in 2:m
        β = _arnoldi_step!(A, V, H, i)
        if β < tol && i > k
            m = i  # happy breakdown
            break
        end
    end

    f = ones(T, m)

    Vₘ = view(V, :, 1:m)
    Hₘ = view(H, 1:m, 1:m)
    qₘ = view(V, :, m + 1)
    βeₘ = view(H, m + 1, 1:m)
    β = real(H[m + 1, m])
    Uₘ = Matrix{T}(I, m, m)
    Uₘᵥ = view(Uₘ, m, 1:m)

    cache0 = similar(v₀, m, m)
    cache1 = similar(v₀, n, m)
    cache2 = similar(H, m)
    sorted_vals = Vector{Int}(undef, m)

    V₁ₖ      = view(V, :, 1:k)
    Vₖ₊₁     = view(V, :, k + 1)
    Hₖ₊₁₁ₖ   = view(H, k + 1, 1:k)
    cache1₁ₖ = view(cache1, :, 1:k)
    cache2₁ₖ = view(cache2, 1:k)

    _update_schur_eigs!(Hₘ, Uₘ, Uₘᵥ, f, k, β, sorted_vals, sortby, rev)

    numops = m
    iter = 0
    while iter < maxiter && count(x -> abs(x) < tol, f) < k && β > tol
        copyto!(cache0, Uₘ)
        mul!(cache1, Vₘ, cache0)
        copyto!(V₁ₖ, cache1₁ₖ)
        copyto!(Vₖ₊₁, qₘ)
        mul!(cache2, transpose(Uₘ), βeₘ)
        copyto!(Hₖ₊₁₁ₖ, cache2₁ₖ)

        @inbounds for j in (k + 1):m
            β = _arnoldi_step!(A, V, H, j)
            if β < tol
                numops += j - k - 1
                break
            end
        end
        _update_schur_eigs!(Hₘ, Uₘ, Uₘᵥ, f, k, β, sorted_vals, sortby, rev)
        numops += m - k - 1
        iter += 1
    end

    Tₘ = Hₘ
    vals = diag(view(Tₘ, 1:k, 1:k))
    VR = _schur_right_eigenvectors(Tₘ, k)
    mul!(cache1₁ₖ, Vₘ, Uₘ * VR)

    idxs = sortperm(vals, by = sortby, rev = rev)
    vals = vals[idxs]
    vecs = cache1₁ₖ[:, idxs]
    converged = iter < maxiter

    return vals, vecs, iter, numops, converged
end

# ---------------------------------------------------------------------
# Concrete eigensolve algorithms
# ---------------------------------------------------------------------

"""
$(TYPEDEF)

Shift-invert Arnoldi with Schur restart, adapted from QuantumToolbox.jl.
No external eigensolver dependency: uses LinearSolve.jl internally for
the inner `(A - σI) x = b` solve.

# Fields
$(TYPEDFIELDS)
"""
struct ArnoldiSI{A, K <: NamedTuple, T} <: EigensolveAlgorithm
    "Number of eigenpairs to extract."
    nev::Int
    "Krylov subspace dimension."
    krylovdim::Int
    "Arnoldi residual tolerance."
    tol::T
    "Maximum number of Arnoldi restarts."
    maxiter::Int
    "LinearSolve.jl algorithm for the inner `(A - σI) x = b` solve."
    linsolve::A
    "Options forwarded to `LinearSolve.init` / `solve!`."
    linsolve_kwargs::K
end

ArnoldiSI(; nev::Int = 1, krylovdim::Int = 20, tol::Real = 1.0e-10,
        maxiter::Int = 200, linsolve = nothing,
        linsolve_kwargs::NamedTuple = (;)) =
    ArnoldiSI(nev, krylovdim, tol, maxiter, linsolve, linsolve_kwargs)

"""
$(TYPEDEF)

Shift-invert inverse-power iteration with Rayleigh-quotient eigenvalue
estimate. Cheapest backend; converges geometrically with ratio
`|λ_target - σ| / |λ_next - σ|`. Slow on clustered spectra.

# Fields
$(TYPEDFIELDS)
"""
struct ShiftInvertPower{A, K <: NamedTuple, T} <: EigensolveAlgorithm
    "Convergence tolerance on the Rayleigh quotient."
    tol::T
    "Maximum number of power-iteration steps."
    maxiter::Int
    "LinearSolve.jl algorithm for the inner `(A - σI) x = b` solve."
    linsolve::A
    "Options forwarded to `LinearSolve.init` / `solve!`."
    linsolve_kwargs::K
end

ShiftInvertPower(; tol::Real = 1.0e-10, maxiter::Int = 1000,
        linsolve = nothing, linsolve_kwargs::NamedTuple = (;)) =
    ShiftInvertPower(tol, maxiter, linsolve, linsolve_kwargs)

"""
$(TYPEDEF)

`KrylovKit.eigsolve`-backed Arnoldi with shift-invert. The
implementation lives in the `CriticalTransitionsKrylovKitExt` package
extension; load `KrylovKit` to enable.

# Fields
$(TYPEDFIELDS)
"""
struct KrylovKitArnoldi{A, K <: NamedTuple, T} <: EigensolveAlgorithm
    "Number of eigenpairs to extract."
    nev::Int
    "Krylov subspace dimension."
    krylovdim::Int
    "Eigsolve tolerance."
    tol::T
    "Maximum number of iterations."
    maxiter::Int
    "LinearSolve.jl algorithm for the inner `(A - σI) x = b` solve."
    linsolve::A
    "Options forwarded to `LinearSolve.init` / `solve!`."
    linsolve_kwargs::K
end

KrylovKitArnoldi(; nev::Int = 1, krylovdim::Int = 30, tol::Real = 1.0e-10,
        maxiter::Int = 200, linsolve = nothing,
        linsolve_kwargs::NamedTuple = (;)) =
    KrylovKitArnoldi(nev, krylovdim, tol, maxiter, linsolve, linsolve_kwargs)

# ---------------------------------------------------------------------
# principal_eigenpair: dispatch on the algorithm struct.
# ---------------------------------------------------------------------

"""
$(TYPEDSIGNATURES)

Compute the eigenpair `(λ, v)` of `A` with eigenvalue closest to `σ`.

`alg` selects the backend ([`ArnoldiSI`](@ref),
[`ShiftInvertPower`](@ref), or `KrylovKitArnoldi`). Each carries its
own LinearSolve.jl `linsolve` algorithm and `linsolve_kwargs` for the
inner `(A - σI) x = b` solve.

Returns `(λ::Complex, v::Vector{Float64}, info::NamedTuple)`. `info`
holds `iter`, `numops`, and `converged`.
"""
function principal_eigenpair end

function principal_eigenpair(A::AbstractMatrix, alg::ArnoldiSI;
        σ::Number = 0.0, v0::Union{Nothing, AbstractVector} = nothing)
    n = size(A, 1)
    T = promote_type(eltype(A), typeof(σ), Float64)
    v_init = v0 === nothing ? normalize!(rand(T, n)) : convert(Vector{T}, v0)
    # Tiny shift away from exact eigenvalues so the inner solve is regular.
    σ_eff = iszero(σ) ? T(-1.0e-12) : σ
    Amap = _build_shift_invert(A, σ_eff, alg.linsolve, alg.linsolve_kwargs)
    k = alg.nev
    m = max(alg.krylovdim, 2k + 1)
    m = min(m, n - 1)
    m < k + 1 && throw(ArgumentError(
        "matrix too small: need n ≥ krylovdim + 1 (got n = $n, krylovdim = $m)",
    ))
    vals_μ, vecs, iter, numops, converged = _eigsolve_arnoldi(
        Amap, v_init, k, m; tol = alg.tol, maxiter = alg.maxiter,
        sortby = abs, rev = true,
    )
    vals_λ = @. σ_eff + 1 / vals_μ
    info = (; iter, numops, converged)
    return vals_λ[1], real.(vecs[:, 1]), info
end

function principal_eigenpair(A::AbstractMatrix, alg::ShiftInvertPower;
        σ::Number = 0.0, v0::Union{Nothing, AbstractVector} = nothing)
    n = size(A, 1)
    T = promote_type(eltype(A), typeof(σ), Float64)
    v = v0 === nothing ? normalize!(rand(T, n)) : copy(convert(Vector{T}, v0))
    normalize!(v)
    σ_eff = iszero(σ) ? T(-1.0e-12) : σ
    Amap = _build_shift_invert(A, σ_eff, alg.linsolve, alg.linsolve_kwargs)
    Av = similar(v)
    λ_prev = T(NaN)
    λ = T(NaN)
    iter = 0
    converged = false
    for outer iter in 1:alg.maxiter
        mul!(Av, Amap, v)
        normalize!(Av)
        Av_raw = A * Av
        λ = dot(Av, Av_raw)
        if !isnan(λ_prev) && abs(λ - λ_prev) <= alg.tol * max(one(real(T)), abs(λ))
            converged = true
            v = Av
            break
        end
        λ_prev = λ
        v = Av
    end
    info = (; iter, numops = iter, converged)
    return Complex(λ), real.(v), info
end
