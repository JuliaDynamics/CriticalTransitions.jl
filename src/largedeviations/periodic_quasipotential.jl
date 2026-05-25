# Periodic Riccati iteration for the quadratic quasi-potential along a limit cycle.

using FastInterpolations: cardinal_interp, PeriodicBC
using LinearAlgebra: Symmetric, mul!, isposdef
using OrdinaryDiffEqLowOrderRK: BS5
using SciMLBase: ODEProblem, solve

"""
    _PRDECache{T}

Caches the τ-grid and `PeriodicBC` reused by every RHS evaluation inside `local_quasipotential`.
Avoids per-step `collect(range(...))` and `PeriodicBC` reconstruction.
"""
struct _PRDECache{T}
    τ_grid::Vector{T}
    bc::PeriodicBC{:exclusive, T, true}
end

function _PRDECache(period::T, Nτ::Int) where {T}
    τ_grid = collect(range(zero(T), period; length = Nτ + 1))[1:Nτ]
    bc = PeriodicBC(endpoint = :exclusive, period = period)
    return _PRDECache{T}(τ_grid, bc)
end

"""
    _prde_rhs!(dG, G, lc::LimitCycleFrame, τ, cache::_PRDECache)

In-place RHS of the periodic Riccati ODE
``\\dot G = -\\tilde M(τ)^\\mathsf{T} G - G \\tilde M(τ) - G \\tilde A(τ) G``.
"""
function _prde_rhs!(dG, G, lc::LimitCycleFrame, τ, cache::_PRDECache)
    period = lc.period
    τw = mod(τ, period)
    M = _interp_matrix(lc.M̃, τw, cache)
    A = _interp_matrix(lc.Ã, τw, cache)
    GA = G * A
    GAG = GA * G
    mul!(dG, transpose(M), G, -1.0, 0.0)
    mul!(dG, G, M, -1.0, 1.0)
    dG .-= GAG
    return dG
end

# Fallback signature without cache (used by single-shot tests / RHS sanity checks).
function _prde_rhs!(dG, G, lc::LimitCycleFrame, τ)
    cache = _PRDECache(lc.period, length(lc.γ))
    return _prde_rhs!(dG, G, lc, τ, cache)
end

# Periodic cubic interpolation of an `m × m × Nτ` array on the precomputed grid.
function _interp_matrix(A3::Array{T, 3}, τw, cache::_PRDECache{T}) where {T}
    m, _, _ = size(A3)
    out = Matrix{T}(undef, m, m)
    for i in 1:m, j in 1:m
        vals = cardinal_interp(cache.τ_grid, view(A3, i, j, :), [τw]; bc = cache.bc)
        out[i, j] = vals[1]
    end
    return out
end

"""
    local_quasipotential(lc::LimitCycleFrame; G0_scale=1e3, maxiters=200, tol=1e-8, alg=BS5())

Solve the periodic Riccati equation `Ġ = -M̃ᵀG - GM̃ - GÃG` for its T-periodic
positive-definite solution `G(τ)`. Returns `G::Array{T, 3}` of shape `(d-1) × (d-1) × Nτ`.
"""
function local_quasipotential(lc::LimitCycleFrame{D, Tf};
        G0_scale::Real = 1.0e3, maxiters::Int = 200, tol::Real = 1.0e-8, alg = BS5(),
    ) where {D, Tf}
    m = D - 1
    Nτ = length(lc.γ)
    period = lc.period
    G_curr = Matrix(Tf(G0_scale) * I, m, m)
    cache = _PRDECache(Tf(period), Nτ)
    f! = (du, u, p, t) -> _prde_rhs!(du, u, lc, t, cache)
    for _ in 1:maxiters
        prob = ODEProblem(f!, G_curr, (zero(Tf), Tf(period)))
        sol = solve(prob, alg; saveat = cache.τ_grid, reltol = 1.0e-9, abstol = 1.0e-12)
        G_next = sol.u[end]
        rel = norm(G_next - G_curr) / max(norm(G_curr), eps(Tf))
        if rel < tol
            G_out = zeros(Tf, m, m, Nτ)
            for k in 1:Nτ
                Gk = sol.u[k]
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
