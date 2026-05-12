module CriticalTransitionsKrylovKitExt

# `KrylovKit.eigsolve` backend for `principal_eigenpair`. Activated when
# the user does `using KrylovKit`.

using CriticalTransitions: KrylovKitArnoldi
import CriticalTransitions: principal_eigenpair
using LinearAlgebra: mul!, normalize!, I
using SparseArrays: sparse
using LinearSolve: LinearProblem, init, solve!
using KrylovKit: eigsolve

# Same `ShiftInvertMap` interface as the in-tree backend, but built
# locally so the extension is self-contained.
struct _SIMap{T, Cache}
    n::Int
    cache::Cache
    buf::Vector{T}
end

Base.eltype(::_SIMap{T}) where {T} = T
Base.size(M::_SIMap) = (M.n, M.n)
Base.size(M::_SIMap, ::Integer) = M.n

function (M::_SIMap)(x::AbstractVector)
    @inbounds @simd for i in 1:M.n
        M.buf[i] = x[i]
    end
    M.cache.b = M.buf
    sol = solve!(M.cache)
    return copy(sol.u)
end

function _build_si(A, σ, linsolve, linsolve_kwargs::NamedTuple)
    n = size(A, 1)
    T = promote_type(eltype(A), typeof(σ), Float64)
    Aσ = sparse(A - σ * I)
    rhs = zeros(T, n)
    prob = LinearProblem{true}(Aσ, rhs)
    cache = linsolve === nothing ?
        init(prob; alias_A = false, alias_b = false, linsolve_kwargs...) :
        init(prob, linsolve; alias_A = false, alias_b = false, linsolve_kwargs...)
    buf = zeros(T, n)
    return _SIMap{T, typeof(cache)}(n, cache, buf)
end

function principal_eigenpair(A::AbstractMatrix, alg::KrylovKitArnoldi;
        σ::Number = 0.0, v0::Union{Nothing, AbstractVector} = nothing)
    n = size(A, 1)
    T = promote_type(eltype(A), typeof(σ), Float64)
    σ_eff = iszero(σ) ? T(-1.0e-12) : σ
    Amap = _build_si(A, σ_eff, alg.linsolve, alg.linsolve_kwargs)
    v_init = v0 === nothing ? normalize!(rand(T, n)) : convert(Vector{T}, v0)
    # `eigsolve` on a function-like map: we want the eigenvalue of A
    # closest to σ, i.e. the largest |μ| of the shift-invert map.
    vals_μ, vecs_μ, info = eigsolve(Amap, v_init, alg.nev, :LM;
        tol = alg.tol, krylovdim = alg.krylovdim, maxiter = alg.maxiter)
    vals_λ = σ_eff .+ 1 ./ vals_μ
    info_nt = (; iter = info.numiter, numops = info.numops, converged = info.converged > 0)
    return vals_λ[1], real.(vecs_μ[1]), info_nt
end

end # module
