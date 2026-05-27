# gMAM-LQA: minimum-action path from a stable limit cycle to a target point.
# Lin, Yu & Zhou (2018), J. Nonlinear Sci. 29:961, Sec. 5.4 eq. (42):
#
#     V(x) = min { Ŝ[φ] + ½ h² ẑᵀ G(τ) ẑ }
#             s.t.  φ(0) = γ(τ) + h Ẽ(τ) ẑ,  φ(1) = x,  |ẑ| = 1
#
# inner: gMAM on the bulk path φ with launch state (τ, ẑ) frozen;
# outer: continuous-τ projected gradient on (τ, ẑ) ∈ S¹ × S^{d-2};
# init:  global scan over Nτ launch-grid candidates (the outer step is local
#        and cannot tunnel between basins of `S(τ)`).

struct _BoundaryWorkspace{Tf}
    path_try::Matrix{Tf}
    alpha::Vector{Tf}
    arc::StepRangeLen{Tf}
    scratch::Vector{Tf}
end
function _BoundaryWorkspace(path::Matrix{Tf}) where {Tf}
    Nt = size(path, 2)
    return _BoundaryWorkspace{Tf}(
        similar(path), zeros(Tf, Nt),
        range(zero(Tf), one(Tf); length = Nt), zeros(Tf, Nt),
    )
end

struct _InnerSweepWS{Tf, R}
    alpha::Vector{Tf}
    arc::R
    scratch::Vector{Tf}
    path_start::Matrix{Tf}
    path_result::Matrix{Tf}
end
function _InnerSweepWS(path::Matrix{Tf}) where {Tf}
    Nt = size(path, 2)
    arc = range(zero(Tf), one(Tf); length = Nt)
    return _InnerSweepWS{Tf, typeof(arc)}(
        zeros(Tf, Nt), arc, zeros(Tf, Nt), similar(path), similar(path),
    )
end

# ---- continuous-τ interpolation ---------------------------------------------
# γ, Ẽ, G are sampled at τ_k = (k-1)·T/Nτ. Linear interpolation lifts them to
# off-grid τ. Across the τ=T seam, anti-periodic columns of Ẽ (where lc.F[j]=-1)
# and entries of G (where lc.S[i,j]=-1) get a sign flip on the k2 side.

@inline function _periodic_lerp_indices(τ::Tf, T_period::Tf, Nτ::Int) where {Tf}
    f = (mod(τ, T_period) / T_period) * Nτ
    k1z = floor(Int, f); α = Tf(f - k1z)
    k1 = mod1(k1z + 1, Nτ); k2 = mod1(k1 + 1, Nτ)
    return k1, k2, α, (k1 == Nτ) && (k2 == 1)
end

@inline function _γ_at(lc::LimitCycleFrame{D, Tf}, τ::Tf) where {D, Tf}
    k1, k2, α, _ = _periodic_lerp_indices(τ, lc.period, length(lc.γ))
    return SVector{D, Tf}(ntuple(i -> (1 - α) * lc.γ[k1][i] + α * lc.γ[k2][i], Val(D)))
end

@inline function _Ẽ_at(lc::LimitCycleFrame{D, Tf, M}, τ::Tf) where {D, Tf, M}
    k1, k2, α, crosses = _periodic_lerp_indices(τ, lc.period, length(lc.γ))
    Ẽ1 = view(lc.Ẽ, :, :, k1); Ẽ2 = view(lc.Ẽ, :, :, k2)
    out = Matrix{Tf}(undef, D, M)
    @inbounds for j in 1:M
        s = (crosses && lc.F[j] == -1) ? -one(Tf) : one(Tf)
        for i in 1:D
            out[i, j] = (1 - α) * Ẽ1[i, j] + α * s * Ẽ2[i, j]
        end
    end
    return out
end

@inline function _G_at(
        G::Array{Tf, 3}, τ::Tf, T_period::Tf, S::SMatrix{M, M, Int8},
    ) where {Tf, M}
    k1, k2, α, crosses = _periodic_lerp_indices(τ, T_period, size(G, 3))
    G1 = view(G, :, :, k1); G2 = view(G, :, :, k2)
    out = Matrix{Tf}(undef, M, M)
    @inbounds for j in 1:M, i in 1:M
        s = (crosses && S[i, j] == -1) ? -one(Tf) : one(Tf)
        out[i, j] = (1 - α) * G1[i, j] + α * s * G2[i, j]
    end
    @inbounds for j in 1:M, i in (j + 1):M
        avg = (out[i, j] + out[j, i]) / 2
        out[i, j] = avg; out[j, i] = avg
    end
    return Symmetric(out)
end

# ---- path primitives --------------------------------------------------------

# Straight line from γ(k0)+h·Ẽ(k0)·ẑ_min to x_f, with ẑ_min = smallest-eigval
# eigenvector of G(k0), signed to bring the tube point closer to x_f.
function _gmam_lqa_straight_init!(
        path::AbstractMatrix{Tf}, lc::LimitCycleFrame{D, Tf, M},
        G::Array{Tf, 3}, x_f, k0::Int, h::Tf,
    ) where {D, Tf, M}
    γ_k0 = collect(lc.γ[k0]); x_f_v = collect(x_f)
    ẑ = eigen(Symmetric(view(G, :, :, k0))).vectors[:, 1]
    Δx = h * (lc.Ẽ[:, :, k0] * ẑ)
    x1 = norm((γ_k0 + Δx) - x_f_v) ≤ norm((γ_k0 - Δx) - x_f_v) ? γ_k0 + Δx : γ_k0 - Δx
    npoints = size(path, 2)
    @inbounds for k in 1:npoints, i in 1:D
        s = (k - 1) / (npoints - 1)
        path[i, k] = (1 - s) * x1[i] + s * x_f_v[i]
    end
    return path
end

function _gmam_lqa_straight_init(
        lc::LimitCycleFrame{D, Tf}, G::Array{Tf, 3}, x_f, k0::Int, h::Tf, npoints::Int,
    ) where {D, Tf}
    return _gmam_lqa_straight_init!(zeros(Tf, D, npoints), lc, G, x_f, k0, h)
end

# Rebuild `dst` from `src` with first node = γ(τ)+h·Ẽ(τ)·ẑ, restore equi-arclength.
@inline function _rebuild_with_launch!(
        dst, src, lc::LimitCycleFrame{D, Tf, M}, τ::Tf, ẑ, h::Tf,
        bws::_BoundaryWorkspace{Tf},
    ) where {D, Tf, M}
    Ẽτ = _Ẽ_at(lc, τ); γτ = _γ_at(lc, τ)
    copyto!(dst, src)
    @inbounds for i in 1:D
        acc = zero(Tf)
        for j in 1:M
            acc += Ẽτ[i, j] * ẑ[j]
        end
        dst[i, 1] = γτ[i] + h * acc
    end
    interpolate_path!(dst, bws.alpha, bws.arc, bws.scratch)
    return dst
end

# ---- inner gMAM sweep -------------------------------------------------------

function _make_path_stepper(path, ws, sys, sw::_InnerSweepWS)
    return function (ϵ)
        geometric_gradient_step!(ws, sys, path; stepsize = ϵ)
        path .= ws.update
        interpolate_path!(path, sw.alpha, sw.arc, sw.scratch)
        return geometric_action(sys, path)
    end
end

function _gmam_lqa_inner_sweep!(
        path, ws, sys, sw::_InnerSweepWS, optimizer::GeometricGradient;
        maxiters::Int = 50,
    )
    step! = _make_path_stepper(path, ws, sys, sw)
    backtracking_optimize!(
        optimizer, step!,
        () -> copyto!(sw.path_start, path), () -> copyto!(path, sw.path_start),
        geometric_action(sys, path); maxiters = maxiters,
    )
    return path
end

function _gmam_lqa_inner_sweep!(
        path, ws, sys, sw::_InnerSweepWS, optimizer::AdaptiveGeometricGradient;
        maxiters::Int = 50,
    )
    step! = _make_path_stepper(path, ws, sys, sw)
    adaptive_optimize!(
        optimizer, step!,
        () -> copyto!(sw.path_start, path), () -> copyto!(path, sw.path_start),
        () -> copyto!(sw.path_result, path), () -> copyto!(path, sw.path_result),
        geometric_action(sys, path); maxiters = maxiters,
    )
    return path
end

# ---- init scan: pick the launch basin ---------------------------------------
# Relax each candidate by `scan_maxiters` inner-sweep iterations and return the
# argmin of the relaxed geometric action. The straight-line action at the tube
# boundary is not a reliable basin indicator (it can be unimodal while the
# relaxed action is multi-modal).

function _pick_initial_k0(
        lc::LimitCycleFrame{D, Tf}, G::Array{Tf, 3}, x_f, sys, optimizer::GMAMOptimizer;
        tube_radius::Tf, npoints::Int, scan_maxiters::Int,
    ) where {D, Tf}
    Nτ = length(lc.γ)
    # maxthreadid (not nthreads) so the cache covers interactive-pool threads too.
    nth = Threads.maxthreadid()
    proto = zeros(Tf, D, npoints)
    caches = [
        (similar(proto), geometric_gradient_workspace(sys, proto), _InnerSweepWS(proto))
            for _ in 1:nth
    ]
    best_S = fill(Tf(Inf), nth)
    best_k = fill(1, nth)
    Threads.@threads for k in 1:Nτ
        tid = Threads.threadid()
        path_tid, ws_tid, sw_tid = caches[tid]
        _gmam_lqa_straight_init!(path_tid, lc, G, x_f, k, tube_radius)
        _gmam_lqa_inner_sweep!(
            path_tid, ws_tid, sys, sw_tid, optimizer; maxiters = scan_maxiters,
        )
        S_k = geometric_action(sys, path_tid)
        if S_k < best_S[tid]
            best_S[tid] = S_k
            best_k[tid] = k
        end
    end
    return best_k[argmin(best_S)]
end

# ---- outer boundary step ----------------------------------------------------

@inline _augmented_F(lc, G, sys, path_consistent, τ::Tf, ẑ, h::Tf) where {Tf} =
    Tf(0.5) * h^2 * dot(ẑ, _G_at(G, τ, lc.period, lc.S) * ẑ) +
    geometric_action(sys, path_consistent)

"""
Continuous-τ projected-gradient step on `(τ, ẑ) ∈ S¹ × S^{d-2}` for
`F(τ, ẑ) = ½ h² ẑᵀ G(τ) ẑ + Ŝ[φ̃(τ, ẑ)]` (Lin et al. 2018 eq. 42). Gradient is
analytic for the LQA term and central-FD for the path term, projected to the
sphere tangent. Backtracking line search at `α ∈ (1, 0.5, 0.1, 0.01)` multipliers
of `(stepsize_τ, stepsize_z)`; defaults make `α=1` shift τ by one grid spacing
`T/Nτ` and ẑ by one sphere quadrant. If no candidate beats the baseline, the
state and path are returned unchanged (monotone by construction).
"""
function _gmam_lqa_boundary_step_continuous(
        state, lc::LimitCycleFrame{D, Tf, M}, G::Array{Tf, 3}, sys, path, x_f, h::Tf,
        bws::_BoundaryWorkspace{Tf};
        stepsize_τ::Tf = lc.period / length(lc.γ), stepsize_z::Tf = one(Tf),
    ) where {D, Tf, M}
    τ_0, ẑ_0 = state
    T = lc.period
    δτ = sqrt(eps(Tf)) * T
    δz = sqrt(eps(Tf))

    # Baseline F (also leaves bws.path_try as the state-consistent path).
    _rebuild_with_launch!(bws.path_try, path, lc, τ_0, ẑ_0, h, bws)
    F_0 = _augmented_F(lc, G, sys, bws.path_try, τ_0, ẑ_0, h)

    # ∂F/∂τ via central FD.
    _rebuild_with_launch!(bws.path_try, path, lc, τ_0 + δτ, ẑ_0, h, bws)
    F_p = _augmented_F(lc, G, sys, bws.path_try, τ_0 + δτ, ẑ_0, h)
    _rebuild_with_launch!(bws.path_try, path, lc, τ_0 - δτ, ẑ_0, h, bws)
    F_m = _augmented_F(lc, G, sys, bws.path_try, τ_0 - δτ, ẑ_0, h)
    grad_τ = (F_p - F_m) / (2δτ)

    # ∂F/∂ẑ: analytic LQA + per-component FD on the path term, projected to
    # the sphere tangent at ẑ_0.
    Gẑ = _G_at(G, τ_0, T, lc.S) * ẑ_0
    grad_z = Vector{Tf}(undef, M)
    ẑ_p = similar(ẑ_0); ẑ_m = similar(ẑ_0)
    @inbounds for i in 1:M
        copyto!(ẑ_p, ẑ_0); ẑ_p[i] += δz
        copyto!(ẑ_m, ẑ_0); ẑ_m[i] -= δz
        _rebuild_with_launch!(bws.path_try, path, lc, τ_0, ẑ_p, h, bws)
        Sp = geometric_action(sys, bws.path_try)
        _rebuild_with_launch!(bws.path_try, path, lc, τ_0, ẑ_m, h, bws)
        Sm = geometric_action(sys, bws.path_try)
        grad_z[i] = h^2 * Gẑ[i] + (Sp - Sm) / (2δz)
    end
    proj = dot(ẑ_0, grad_z)
    @inbounds for i in 1:M
        grad_z[i] -= proj * ẑ_0[i]
    end

    # Normalize so step α = 1 has magnitude (stepsize_τ, stepsize_z).
    dir_τ = grad_τ / max(abs(grad_τ), eps(Tf))
    dir_z = grad_z ./ max(norm(grad_z), eps(Tf))

    # Backtracking line search.
    best_state = (mod(τ_0, T), ẑ_0); best_F = F_0; best_path = nothing
    for α in (one(Tf), Tf(0.5), Tf(0.1), Tf(0.01))
        τ_try = τ_0 - α * stepsize_τ * dir_τ
        ẑ_try = ẑ_0 .- α * stepsize_z .* dir_z
        nrm = norm(ẑ_try); nrm < eps(Tf) && continue
        ẑ_try ./= nrm
        _rebuild_with_launch!(bws.path_try, path, lc, τ_try, ẑ_try, h, bws)
        F_try = _augmented_F(lc, G, sys, bws.path_try, τ_try, ẑ_try, h)
        if F_try < best_F
            best_state = (mod(τ_try, T), ẑ_try); best_F = F_try
            best_path = copy(bws.path_try)
        end
    end
    return best_state, best_path === nothing ? path : best_path
end

# Convenience form (used by tests): build the workspace per call.
_gmam_lqa_boundary_step_continuous(state, lc, G, sys, path, x_f, h) =
    _gmam_lqa_boundary_step_continuous(
    state, lc, G, sys, path, x_f, eltype(path)(h), _BoundaryWorkspace(path),
)

# ---- public entry -----------------------------------------------------------

"""
    minimize_geometric_action(sys::CoupledSDEs, lc::LimitCycleFrame, x_f,
        optimizer::GMAMOptimizer = GeometricGradient(; stepsize = 0.01); kwargs...)

gMAM-LQA: minimum-action path from a stable limit cycle `Γ` (encoded by `lc`) to
`x_f`, after Lin, Yu & Zhou (2018), J. Nonlinear Sci. 29:961, Sec. 5.4 eq. (42).
Inside a tube of radius `h` around `Γ` the quasi-potential is approximated by the
quadratic form `½ ẑᵀ G(τ) ẑ` (from [`local_quasipotential`](@ref)); outside the
tube `Ŝ[φ]` is minimized by gMAM with the launch end pinned to the tube surface.
The returned action is `Ŝ[φ] + ½ h² ẑᵀ G(τ_start) ẑ`, with `(τ_start, ẑ)` jointly
optimized on `S¹ × S^{d-2}` (continuous τ).

Pass `GeometricGradient` for backtracking step-size control or
`AdaptiveGeometricGradient` for probe-based adaptation (more robust on
underdamped systems).

## Keyword arguments

  - `G::Array{T, 3} = local_quasipotential(lc)`: transverse LQA Hessian along `Γ`.
  - `tube_radius::Real = 0.05`: tube radius `h`.
  - `npoints::Int = 100`: path discretization size.
  - `maxiters::Int = 50`: outer iterations.
  - `inner_maxiters::Int = 50`: inner-sweep iterations per outer iter.
  - `abstol::Real = NaN`, `reltol::Real = 1e-4`: outer convergence tolerances.
  - `k0_init::Union{Nothing, Int} = nothing`: initial launch grid index. Default
    `nothing` scans all `Nτ` candidates with `scan_maxiters` inner-sweep
    iterations each and picks the basin (essential for multi-modal `S(τ)`). Pass
    an integer to skip the scan.
  - `scan_maxiters::Int = inner_maxiters`: per-candidate iterations during the
    init scan. Lower (e.g. `100`) for cheaper init on unimodal landscapes.
  - `verbose::Bool = false`, `show_progress::Bool = false`.
"""
function minimize_geometric_action(
        sys::CoupledSDEs, lc::LimitCycleFrame{D, Tf}, x_f,
        optimizer::GMAMOptimizer = GeometricGradient(; stepsize = 0.01);
        G::Array{Tf, 3} = local_quasipotential(lc),
        tube_radius::Real = 0.05, npoints::Int = 100,
        maxiters::Int = 50, inner_maxiters::Int = 50,
        abstol::Real = NaN, reltol::Real = 1.0e-4,
        verbose::Bool = false, show_progress::Bool = false,
        k0_init::Union{Nothing, Int} = nothing, scan_maxiters::Int = inner_maxiters,
    ) where {D, Tf}
    h = Tf(tube_radius)

    # Init: global scan over launch candidates (or user-supplied k0).
    k0 = k0_init === nothing ?
        _pick_initial_k0(
            lc, G, x_f, sys, optimizer;
            tube_radius = h, npoints = npoints, scan_maxiters = scan_maxiters
        ) :
        k0_init
    path = _gmam_lqa_straight_init(lc, G, x_f, k0, h, npoints)

    # Workspaces + continuous-τ launch state recovered from grid index k0.
    ws = geometric_gradient_workspace(sys, path)
    bws = _BoundaryWorkspace(path)
    sw = _InnerSweepWS(path)
    z_raw = lc.Ẽ[:, :, k0]' * (path[:, 1] - collect(lc.γ[k0]))
    state = (Tf((k0 - 1) * lc.period / length(lc.γ)), z_raw / norm(z_raw))

    # Outer loop.
    S = _augmented_F(lc, G, sys, path, state[1], state[2], h)
    progress = Progress(maxiters; dt = 0.5, enabled = show_progress)
    for it in 1:maxiters
        _gmam_lqa_inner_sweep!(path, ws, sys, sw, optimizer; maxiters = inner_maxiters)
        state, path = _gmam_lqa_boundary_step_continuous(state, lc, G, sys, path, x_f, h, bws)
        S_new = _augmented_F(lc, G, sys, path, state[1], state[2], h)
        abs_change = abs(S_new - S)
        rel_change = abs_change / max(abs(S_new), eps(Tf))
        verbose && @info "gMAM-LQA iter $(it): S=$(S_new), abs=$(abs_change), rel=$(rel_change)"
        S = S_new
        if (isfinite(abstol) && abs_change < abstol) ||
                (isfinite(reltol) && rel_change < reltol)
            break
        end
        next!(
            progress; showvalues = [
                ("iter", it), ("S", round(S; sigdigits = 6)),
                ("rel", round(rel_change; sigdigits = 3)),
            ]
        )
    end
    return MinimumActionPath(StateSpaceSet(Matrix(path')), S)
end
