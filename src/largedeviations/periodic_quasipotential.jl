# Periodic Riccati iteration for the quadratic quasi-potential along a limit cycle.
#
# Supports both purely-periodic and antiperiodic transverse Floquet bundles. In the
# antiperiodic case (any `lc.F[j] == -1`) the period map from `G_n` to `G_{n+1} =
# G(2T)` decomposes into two legs that share the same RHS and the same cache, with
# a free elementwise sign-mask transform between them. This avoids any τ-dependent
# branch in the hot RHS path.

using FastInterpolations: cardinal_interp, PeriodicBC
using LinearAlgebra: Symmetric, mul!, isposdef
using OrdinaryDiffEqLowOrderRK: BS5
using SciMLBase: ODEProblem, remake, solve
using StaticArrays: SMatrix

const _SMATRIX_MAX_M = 4

struct _PRDECache{T, ITP}
    τ_grid::Vector{T}
    M_itps::Matrix{ITP}
    A_itps::Matrix{ITP}
    M_buf::Matrix{T}
    A_buf::Matrix{T}
    GA_buf::Matrix{T}
end

function _PRDECache(lc::LimitCycleFrame{D, T}) where {D, T}
    Nτ = length(lc.γ)
    m = D - 1
    period = lc.period
    τ_grid = collect(range(zero(T), period; length = Nτ + 1))[1:Nτ]
    bc = PeriodicBC(endpoint = :exclusive, period = period)
    sample_itp = cardinal_interp(τ_grid, view(lc.M̃, 1, 1, :); bc)
    ITP = typeof(sample_itp)
    M_itps = Matrix{ITP}(undef, m, m)
    A_itps = Matrix{ITP}(undef, m, m)
    for i in 1:m, j in 1:m
        M_itps[i, j] = cardinal_interp(τ_grid, view(lc.M̃, i, j, :); bc)
        A_itps[i, j] = cardinal_interp(τ_grid, view(lc.Ã, i, j, :); bc)
    end
    return _PRDECache{T, ITP}(
        τ_grid, M_itps, A_itps,
        Matrix{T}(undef, m, m), Matrix{T}(undef, m, m), Matrix{T}(undef, m, m),
    )
end

@inline function _eval_itps!(out::Matrix{T}, itps::Matrix, τw::T) where {T}
    m = size(out, 1)
    @inbounds for i in 1:m, j in 1:m
        out[i, j] = itps[i, j](τw)
    end
    return out
end

# Matrix-path RHS: Ġ = -M'·G - G·M - G·Ã·G. Symmetry of G is preserved analytically
# and re-imposed at the outer-iteration boundary; no per-step projection needed.
function _prde_rhs!(dG, G, lc::LimitCycleFrame, τ, cache::_PRDECache)
    τw = mod(τ, lc.period)
    M = _eval_itps!(cache.M_buf, cache.M_itps, τw)
    A = _eval_itps!(cache.A_buf, cache.A_itps, τw)
    mul!(cache.GA_buf, G, A)
    mul!(dG, transpose(M), G, -1, 0)
    mul!(dG, G, M, -1, 1)
    mul!(dG, cache.GA_buf, G, -1, 1)
    return dG
end

function _prde_rhs!(dG, G, lc::LimitCycleFrame, τ)
    return _prde_rhs!(dG, G, lc, τ, _PRDECache(lc))
end

# SMatrix-path RHS for small m (≤ _SMATRIX_MAX_M). Fully unrolled by the compiler;
# no BLAS dispatch and no heap allocation.
@generated function _eval_itp_smatrix(itps::Matrix, τw::T, ::Val{m}) where {T, m}
    exprs = [:(itps[$i, $j](τw)) for j in 1:m for i in 1:m]
    return :(SMatrix{$m, $m, $T}($(exprs...)))
end

@inline function _prde_rhs_sa(G::SMatrix{m, m, T}, lc, τ, cache::_PRDECache, ::Val{m}) where {m, T}
    τw = mod(τ, lc.period)
    M = _eval_itp_smatrix(cache.M_itps, τw, Val(m))
    A = _eval_itp_smatrix(cache.A_itps, τw, Val(m))
    return -transpose(M) * G - G * M - G * A * G
end

# Apply the elementwise sign-mask `S` (entries ±1) in place: this is the action of
# `diag(F) · X · diag(F)` for free.
@inline function _apply_sign_mask!(X::AbstractMatrix, S::SMatrix)
    @inbounds for j in 1:size(X, 2), i in 1:size(X, 1)
        X[i, j] *= S[i, j]
    end
    return X
end

@inline _apply_sign_mask(X::SMatrix{m, m, T}, S::SMatrix{m, m, Int8}) where {m, T} = T.(S) .* X

"""
    local_quasipotential(lc::LimitCycleFrame; G0=nothing, G0_scale=1e3,
                         maxiters=200, tol=1e-8, alg=BS5())

Solve the periodic Riccati equation `Ġ = -M̃ᵀG - GM̃ - GÃG` for its positive-definite
periodic solution `G(τ)` along the limit cycle `lc`. Returns `G::Array{T, 3}` of shape
`(d-1) × (d-1) × Nτ` storing `G(τ_k)` for `τ_k ∈ [0, T)`.

When `all(lc.F .== 1)` (purely periodic transverse bundle), `G` is `T`-periodic and is
returned directly. When `any(lc.F .== -1)` (antiperiodic case), `G` is `2T`-periodic
with `G(τ + T) = diag(F) · G(τ) · diag(F)`; the second half is recoverable on the fly
via [`G_at`](@ref).

Pass `G0::Union{Nothing, AbstractArray}` to warm-start the outer iteration from a
previous result; for parameter sweeps this typically reduces the iteration count by
roughly an order of magnitude.
"""
function local_quasipotential(
        lc::LimitCycleFrame{D, Tf, M};
        G0 = nothing,
        G0_scale::Real = 1.0e3, maxiters::Int = 200, tol::Real = 1.0e-8, alg = BS5(),
    ) where {D, Tf, M}
    Nτ = length(lc.γ)
    period = lc.period
    cache = _PRDECache(lc)
    G_curr_mat = _init_G(Tf, M, G0, G0_scale)
    antiperiodic = any(lc.F .== -1)
    if M ≤ _SMATRIX_MAX_M
        return _local_quasipotential_sa(
            lc, cache, G_curr_mat, alg, period, Nτ, maxiters, tol, antiperiodic, Val(M)
        )
    else
        return _local_quasipotential_mat(
            lc, cache, G_curr_mat, alg, period, Nτ, maxiters, tol, antiperiodic
        )
    end
end

function _init_G(::Type{Tf}, M::Int, G0, G0_scale::Real) where {Tf}
    if G0 === nothing
        return Matrix(Tf(G0_scale) * I, M, M)
    elseif ndims(G0) == 3
        return Matrix{Tf}(view(G0, :, :, 1))
    else
        return Matrix{Tf}(G0)
    end
end

function _local_quasipotential_sa(
        lc::LimitCycleFrame{D, Tf, M}, cache, G_curr_mat, alg,
        period, Nτ, maxiters, tol, antiperiodic, ::Val{m},
    ) where {D, Tf, M, m}
    G_curr = SMatrix{m, m, Tf}(G_curr_mat)
    S_F = lc.S
    f = (u, p, t) -> _prde_rhs_sa(u, lc, t, cache, Val(m))
    prob = ODEProblem(f, G_curr, (zero(Tf), Tf(period)))
    for _ in 1:maxiters
        sol1 = solve(
            remake(prob; u0 = G_curr), alg;
            saveat = cache.τ_grid, reltol = 1.0e-9, abstol = 1.0e-12,
        )
        G_T = sol1.u[end]
        if antiperiodic
            G̃_0 = _apply_sign_mask(G_T, S_F)
            sol2 = solve(
                remake(prob; u0 = G̃_0), alg;
                reltol = 1.0e-9, abstol = 1.0e-12,
            )
            G_next = _apply_sign_mask(sol2.u[end], S_F)
        else
            G_next = G_T
        end
        rel = norm(G_next - G_curr) / max(norm(G_curr), eps(Tf))
        if rel < tol
            G_out = zeros(Tf, m, m, Nτ)
            for k in 1:Nτ
                Gk = sol1.u[k]
                G_out[:, :, k] .= Matrix((Gk + transpose(Gk)) / 2)
            end
            for k in 1:Nτ
                isposdef(Symmetric(G_out[:, :, k])) || error(
                    "local_quasipotential: returned G is not positive definite at τ index $(k); rel=$(rel).",
                )
            end
            return G_out
        end
        G_curr = (G_next + transpose(G_next)) / 2
    end
    return error("local_quasipotential: failed to converge in $(maxiters) outer iterations.")
end

function _local_quasipotential_mat(
        lc::LimitCycleFrame{D, Tf, M}, cache, G_curr, alg,
        period, Nτ, maxiters, tol, antiperiodic,
    ) where {D, Tf, M}
    S_F = lc.S
    f! = (du, u, p, t) -> _prde_rhs!(du, u, lc, t, cache)
    prob = ODEProblem(f!, G_curr, (zero(Tf), Tf(period)))
    for _ in 1:maxiters
        sol1 = solve(
            remake(prob; u0 = G_curr), alg;
            saveat = cache.τ_grid, reltol = 1.0e-9, abstol = 1.0e-12,
        )
        G_T = copy(sol1.u[end])
        if antiperiodic
            _apply_sign_mask!(G_T, S_F)
            sol2 = solve(
                remake(prob; u0 = G_T), alg;
                reltol = 1.0e-9, abstol = 1.0e-12,
            )
            G_next = copy(sol2.u[end])
            _apply_sign_mask!(G_next, S_F)
        else
            G_next = G_T
        end
        rel = norm(G_next - G_curr) / max(norm(G_curr), eps(Tf))
        if rel < tol
            G_out = zeros(Tf, M, M, Nτ)
            for k in 1:Nτ
                Gk = sol1.u[k]
                G_out[:, :, k] = (Gk + Gk') / 2
            end
            for k in 1:Nτ
                isposdef(Symmetric(G_out[:, :, k])) || error(
                    "local_quasipotential: returned G is not positive definite at τ index $(k); rel=$(rel).",
                )
            end
            return G_out
        end
        G_curr = (G_next + G_next') / 2
    end
    return error("local_quasipotential: failed to converge in $(maxiters) outer iterations.")
end

"""
    G_at(lc::LimitCycleFrame, G::Array, τ::Real)

Return `G(τ)` for any `τ ∈ [0, 2T)`. For `τ ∈ [0, T)` this is a nearest-grid-point
lookup into the stored `G`. For `τ ∈ [T, 2T)` (only meaningful in the antiperiodic
case) the second-half value is generated on the fly via the identity
`G(τ + T) = diag(F) · G(τ) · diag(F)`, applied as a Hadamard product with `lc.S`.
"""
function G_at(lc::LimitCycleFrame{D, T, M}, G::Array{T, 3}, τ::Real) where {D, T, M}
    Nτ = size(G, 3)
    Δτ = lc.period / Nτ
    τw = mod(τ, 2 * lc.period)
    if τw < lc.period
        k = mod1(floor(Int, τw / Δτ) + 1, Nτ)
        return SMatrix{M, M, T}(view(G, :, :, k))
    else
        k = mod1(floor(Int, (τw - lc.period) / Δτ) + 1, Nτ)
        Gk = SMatrix{M, M, T}(view(G, :, :, k))
        return T.(lc.S) .* Gk
    end
end
