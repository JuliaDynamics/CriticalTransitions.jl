# gMAM-LQA: minimum-action path anchored on a limit-cycle tube surface.

using LinearAlgebra: eigen, dot, norm, Symmetric

function _γ_vectors(lc::LimitCycleFrame{D, Tf}) where {D, Tf}
    Nτ = length(lc.γ)
    out = Vector{Vector{Tf}}(undef, Nτ)
    for k in 1:Nτ
        out[k] = collect(lc.γ[k])
    end
    return out
end

function _symmetrize_G(G::Array{Tf, 3}) where {Tf}
    Nτ = size(G, 3)
    return [Symmetric((G[:, :, k] + G[:, :, k]') / 2) for k in 1:Nτ]
end

"""
Build the straight-line initial guess for gMAM-LQA: first node at `γ(τ_start) +
h·Ẽ(τ_start)·ẑ` where `τ_start` is the orbit grid point closest to `x_f` and `ẑ` is
the smallest-eigenvalue eigenvector of `G(τ_start)`. Eigenvectors are determined up
to sign; we pick the sign that places the tube point closer to `x_f`. Last node is `x_f`.
"""
function _gmam_lqa_initial_path(
        lc::LimitCycleFrame, G::Array{Tf, 3}, x_f;
        tube_radius::Real, npoints::Int,
    ) where {Tf}
    γs = _γ_vectors(lc)
    x_f_v = collect(x_f)
    k0 = argmin(k -> norm(γs[k] - x_f_v), eachindex(γs))
    F = eigen(Symmetric((G[:, :, k0] + G[:, :, k0]') / 2))
    ẑ = F.vectors[:, 1]
    h = Tf(tube_radius)
    x1_plus = γs[k0] + h * (lc.Ẽ[:, :, k0] * ẑ)
    x1_minus = γs[k0] - h * (lc.Ẽ[:, :, k0] * ẑ)
    x1 = norm(x1_plus - x_f_v) ≤ norm(x1_minus - x_f_v) ? x1_plus : x1_minus
    d = length(x1)
    path = zeros(Tf, d, npoints)
    for k in 1:npoints
        s = (k - 1) / (npoints - 1)
        path[:, k] = (1 - s) * x1 + s * x_f_v
    end
    return path
end

"""
Preallocated buffers used by `_gmam_lqa_boundary_step` across its 9-candidate inner loop
to avoid per-iteration allocation churn.
"""
struct _BoundaryWorkspace{Tf}
    path_try::Matrix{Tf}
    alpha::Vector{Tf}
    arc::StepRangeLen{Tf}
    scratch::Vector{Tf}
end

function _BoundaryWorkspace(path::Matrix{Tf}) where {Tf}
    Nt = size(path, 2)
    return _BoundaryWorkspace{Tf}(
        similar(path),
        zeros(Tf, Nt),
        range(zero(Tf), one(Tf); length = Nt),
        zeros(Tf, Nt),
    )
end

"""
Projected-gradient step on `(k0, ẑ)` for the augmented objective
`S(k0, ẑ) = ½ h² ẑᵀ G(τ_{k0}) ẑ + geometric_action(sys, path_with_new_first_node)`.
Searches over `k0 ∈ {k0, k0±1}` and along the spherical-tangent gradient for `ẑ`
with backtracking step sizes `(1.0, 0.5, 0.1)`. Uses `bws` for scratch.
"""
function _gmam_lqa_boundary_step(
        state, lc::LimitCycleFrame, G_sym::Vector{<:Symmetric},
        γs::Vector{Vector{Tf}}, sys, path, x_f, h, bws::_BoundaryWorkspace{Tf}
    ) where {Tf}
    Nτ = length(γs)
    k0, ẑ = state
    candidates_k = (k0, mod1(k0 + 1, Nτ), mod1(k0 - 1, Nτ))
    best_S = Tf(Inf)
    best_k = k0
    best_ẑ = ẑ
    best_path = path
    for kc in candidates_k
        Gk = G_sym[kc]
        grad_local = h^2 * (Gk * ẑ)
        g_t = grad_local - dot(ẑ, grad_local) * ẑ
        for step in (one(Tf), Tf(0.5), Tf(0.1))
            ẑ_try = ẑ - step * g_t
            ẑ_try ./= norm(ẑ_try)
            x1_try = γs[kc] + h * (lc.Ẽ[:, :, kc] * ẑ_try)
            copyto!(bws.path_try, path)
            bws.path_try[:, 1] .= x1_try
            interpolate_path!(bws.path_try, bws.alpha, bws.arc, bws.scratch)
            S = Tf(0.5) * h^2 * dot(ẑ_try, Gk * ẑ_try) + geometric_action(sys, bws.path_try)
            if S < best_S
                best_S = S; best_k = kc; best_ẑ = ẑ_try
                best_path = copy(bws.path_try)
            end
        end
    end
    return (best_k, best_ẑ), best_path
end

function _gmam_lqa_boundary_step(state, lc::LimitCycleFrame, G::Array{Tf, 3}, sys, path, x_f, h) where {Tf}
    return _gmam_lqa_boundary_step(state, lc, _symmetrize_G(G), _γ_vectors(lc), sys, path, x_f, h, _BoundaryWorkspace(path))
end

"""
Run `N_inner` `geometric_gradient_step!` updates on `path`, keeping the first and last
columns fixed via arclength reinterpolation. Mutates `path` in place.
"""
function _gmam_lqa_inner_sweep!(path, ws, sys; N_inner::Int = 50, stepsize::Real = 1.0)
    Nt = size(path, 2)
    Tf = eltype(path)
    alpha = zeros(Tf, Nt)
    arc = range(zero(Tf), one(Tf); length = Nt)
    scratch = similar(alpha)
    for _ in 1:N_inner
        geometric_gradient_step!(ws, sys, path; stepsize)
        path .= ws.update
        interpolate_path!(path, alpha, arc, scratch)
    end
    return path
end

"""
    minimize_geometric_action(sys::CoupledSDEs, lc::LimitCycleFrame, x_f;
        G = local_quasipotential(lc), tube_radius = 0.05, npoints = 100,
        N_inner = 50, maxiters = 50, reltol = 1e-4, stepsize = 0.01, verbose = false)

gMAM-LQA: minimum action path from limit cycle Γ (encoded by `lc`) to `x_f`, using the
local quadratic approximation of the quasi-potential inside a tube of radius `tube_radius`
around Γ. Returns a [`MinimumActionPath`](@ref) with the augmented action
`Ŝ[φ] + ½ h² ẑᵀ G(τ_start) ẑ`.

Pass `G` explicitly when calling repeatedly to avoid re-running `local_quasipotential` each call.
"""
function minimize_geometric_action(
        sys::CoupledSDEs, lc::LimitCycleFrame{D, Tf}, x_f;
        G::Array{Tf, 3} = local_quasipotential(lc),
        tube_radius::Real = 0.05,
        npoints::Int = 100,
        N_inner::Int = 50,
        maxiters::Int = 50,
        reltol::Real = 1.0e-4,
        stepsize::Real = 0.01,
        verbose::Bool = false,
    ) where {D, Tf}
    h = Tf(tube_radius)
    path = _gmam_lqa_initial_path(lc, G, x_f; tube_radius = h, npoints = npoints)
    ws = geometric_gradient_workspace(sys, path)
    γs = _γ_vectors(lc)
    G_sym = _symmetrize_G(G)
    bws = _BoundaryWorkspace(path)
    k0 = argmin(k -> norm(γs[k] - path[:, 1]), eachindex(γs))
    z_raw = lc.Ẽ[:, :, k0]' * (path[:, 1] - γs[k0])
    ẑ = z_raw / norm(z_raw)
    state = (k0, ẑ)
    S = Tf(0.5) * h^2 * dot(ẑ, G_sym[k0] * ẑ) + geometric_action(sys, path)
    for it in 1:maxiters
        _gmam_lqa_inner_sweep!(path, ws, sys; N_inner = N_inner, stepsize = stepsize)
        state, path = _gmam_lqa_boundary_step(state, lc, G_sym, γs, sys, path, x_f, h, bws)
        S_new = Tf(0.5) * h^2 * dot(state[2], G_sym[state[1]] * state[2]) + geometric_action(sys, path)
        rel = abs(S_new - S) / max(abs(S_new), eps(Tf))
        verbose && @info "gMAM-LQA iter $(it): S=$(S_new), rel=$(rel)"
        S = S_new
        rel < reltol && break
    end
    return MinimumActionPath(StateSpaceSet(Matrix(path')), S)
end
