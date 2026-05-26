# Periodic Riccati iteration for the quadratic quasi-potential along a limit cycle.

using FastInterpolations: cardinal_interp, PeriodicBC
using LinearAlgebra: Symmetric, mul!, isposdef
using OrdinaryDiffEqLowOrderRK: BS5
using SciMLBase: ODEProblem, remake, solve

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

function _eval_itps!(out::Matrix{T}, itps::Matrix, τw::T) where {T}
    m = size(out, 1)
    @inbounds for i in 1:m, j in 1:m
        out[i, j] = itps[i, j](τw)
    end
    return out
end

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

"""
    local_quasipotential(lc::LimitCycleFrame; G0_scale=1e3, maxiters=200, tol=1e-8, alg=BS5())

Solve the periodic Riccati equation `Ġ = -M̃ᵀG - GM̃ - GÃG` for its T-periodic
positive-definite solution `G(τ)`. Returns `G::Array{T, 3}` of shape `(d-1) × (d-1) × Nτ`.
"""
function local_quasipotential(
        lc::LimitCycleFrame{D, Tf};
        G0_scale::Real = 1.0e3, maxiters::Int = 200, tol::Real = 1.0e-8, alg = BS5(),
    ) where {D, Tf}
    m = D - 1
    Nτ = length(lc.γ)
    period = lc.period
    G_curr = Matrix(Tf(G0_scale) * I, m, m)
    cache = _PRDECache(lc)
    f! = (du, u, p, t) -> _prde_rhs!(du, u, lc, t, cache)
    prob = ODEProblem(f!, G_curr, (zero(Tf), Tf(period)))
    for _ in 1:maxiters
        sol = solve(remake(prob; u0 = G_curr), alg; saveat = cache.τ_grid, reltol = 1.0e-9, abstol = 1.0e-12)
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
