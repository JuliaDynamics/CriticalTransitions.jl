"""
    minimize_geometric_action(sys::CoupledSDEs, x_i, x_f, optimizer=GeometricGradient(); kwargs...)
    minimize_geometric_action(sys::CoupledSDEs, init::AbstractMatrix, optimizer=GeometricGradient(); kwargs...)

Computes the minimizer of the geometric Freidlin-Wentzell action based on the geometric
minimum action method (gMAM), using optimizers of Optimization.jl or the original
formulation by [heymann_pathways_2008](@citet) (projected gradient descent with optional
backtracking).

The minimizer is computed for system `sys` over all paths from `x_i` to `x_f`. To set an
initial path different from a straight line, see the multiple-dispatch method
`minimize_geometric_action(sys::CoupledSDEs, init::AbstractMatrix, optimizer; kwargs...)`.

Returns a [`MinimumActionPath`](@ref).

## Keyword arguments
  - `npoints::Int = 100`: number of discretization points (endpoint form only)
  - `maxiters::Int = 100`: maximum optimizer iterations
  - `abstol::Real = NaN`, `reltol::Real = NaN`: action-change convergence criteria
  - `ad_type = OptimizationBase.AutoFiniteDiff()` (Optimization.jl path only)
  - `verbose::Bool = false`, `show_progress::Bool = true`
"""
function minimize_geometric_action(
        sys::CoupledSDEs, x_i, x_f, optimizer = GeometricGradient();
        npoints::Int = 100, kwargs...,
    )
    path = reduce(hcat, range(x_i, x_f; length = npoints))
    return minimize_geometric_action(sys, path, optimizer; kwargs...)
end

minimize_geometric_action(sys::CoupledSDEs, init::AbstractMatrix; kwargs...) =
    minimize_geometric_action(sys, init, GeometricGradient(); kwargs...)

# Shared setup for both backtracking and Optimization.jl paths.
function _gmam_setup(sys::CoupledSDEs, init)
    proper_FW_system(sys)
    _validate_and_classify_a(_trace_normalized_a(sys), collect(view(init, :, 1)))
    N = size(init, 2)
    T = eltype(init)
    alpha = zeros(T, N)
    arc = range(0, 1.0; length = N)
    A_at = _action_metric(sys)
    b_fn = let sys = sys
        x -> drift(sys, x)
    end
    x_i = init[:, 1]; x_f = init[:, end]
    v_buf = similar(init)
    integrand_buf = zeros(T, N)
    S = let b_fn = b_fn, x_i = x_i, x_f = x_f, A_at = A_at, v_buf = v_buf, integrand_buf = integrand_buf
        x -> _geometric_action_from_drift!(b_fn, fix_ends(x, x_i, x_f), 1.0, A_at, v_buf, integrand_buf)
    end
    return alpha, arc, S
end

function minimize_geometric_action(
        sys::CoupledSDEs, init::AbstractMatrix, optimizer::GeometricGradient;
        maxiters::Int = 100, abstol::Real = NaN, reltol::Real = NaN,
        verbose = false, show_progress = true,
    )
    alpha, arc, S = _gmam_setup(sys, init)
    path = Matrix(init)
    ws = geometric_gradient_workspace(sys, path)
    path_prev = similar(path)
    interp_scratch = similar(alpha)
    interpolate_path!(path, alpha, arc, interp_scratch)
    initial_action = S(path)

    function try_step!(ϵ)
        geometric_gradient_step!(ws, sys, path; stepsize = ϵ)
        path .= ws.update
        interpolate_path!(path, alpha, arc, interp_scratch)
        return S(path)
    end
    save!() = copyto!(path_prev, path)
    restore!() = copyto!(path, path_prev)

    current_action, _ = backtracking_optimize!(
        optimizer, try_step!, save!, restore!, initial_action;
        maxiters, abstol, reltol, verbose, show_progress,
    )
    return MinimumActionPath(StateSpaceSet(path'), current_action)
end

function minimize_geometric_action(
        sys::CoupledSDEs, init::AbstractMatrix, optimizer;
        maxiters::Int = 100, abstol::Real = NaN, reltol::Real = NaN,
        ad_type = OptimizationBase.AutoFiniteDiff(),
        show_progress = true,
    )
    init = Matrix(init)
    alpha, arc, S = _gmam_setup(sys, init)
    interp_scratch = similar(alpha)
    optf = SciMLBase.OptimizationFunction((x, _) -> S(x), ad_type)
    prob = SciMLBase.OptimizationProblem(optf, init, ())
    prog = Progress(maxiters; enabled = show_progress)
    function callback(state, _)
        interpolate_path!(state.u, alpha, arc, interp_scratch)
        show_progress ? next!(prog) : nothing
        return false
    end
    sol = solve(prob, optimizer; maxiters, callback, abstol, reltol)
    return MinimumActionPath(StateSpaceSet(sol.u'), S(sol.u))
end

struct GeometricGradientWorkspace{T, A, F, JC, LC, LU}
    update::Matrix{T}
    x_prime::Matrix{T}
    lambdas::Vector{T}
    lambdas_prime::Vector{T}
    drift_cache::Matrix{T}
    theta_cache::Matrix{T}
    rhs_explicit::Matrix{T}
    arc::StepRangeLen{T}
    a_func::A
    a_const_lu::LU            # LU(a) for constant non-diagonal a; `nothing` otherwise
    jac_fn::F                # x -> f(x, p, 0.0)
    jac_cfg::JC              # ForwardDiff.JacobianConfig (preallocated)
    J_buf::Matrix{T}         # output Jacobian buffer, Nx×Nx
    x_buf::Vector{T}         # input x buffer for jacobian!, Nx
    linear_cache::LC
    fd_probe::Vector{T}     # FD probe direction; size Nx
    Hphi::Vector{T}         # H_φ accumulator; size Nx
    Hpphi::Vector{T}        # H_pφ accumulator; size Nx
    aHphi::Vector{T}        # a · H_φ; size Nx
    Ainv_b::Vector{T}       # a⁻¹·b scratch for _lambda_theta!; size Nx
    Ainv_φp::Vector{T}      # a⁻¹·φ' scratch for _lambda_theta!; size Nx
    x_probe::Vector{T}      # x ± h·eₗ scratch for state-dep ∂a; size Nx
    dla_buf::Matrix{T}      # ∂_l a(x) scratch; size Nx×Nx
end

function geometric_gradient_workspace(sys::CoupledSDEs, path)
    T = eltype(path)
    Nx, N = size(path)
    proper_FW_system(sys)
    a_func = _trace_normalized_a(sys)
    _validate_and_classify_a(a_func, collect(view(path, :, 1)))
    params = current_parameters(sys)
    f_rule = dynamic_rule(sys)
    jac_fn = let f = f_rule, params = params
        x -> f(x, params, 0.0)
    end
    x_buf = collect(view(path, :, 1))
    J_buf = zeros(T, Nx, Nx)
    jac_cfg = ForwardDiff.JacobianConfig(jac_fn, x_buf)
    linear_cache = _init_tridiag_cache(T, N)
    # Precompute LU(a) for constant non-diagonal a. (Diagonal a uses the scalar fast path.)
    a_const_lu = let a0 = a_func(view(path, :, 1))
        if is_constant(a_func) === Val(true) && !(a0 isa LinearAlgebra.Diagonal)
            LinearAlgebra.lu(a0 isa Matrix ? Matrix{T}(a0) : Matrix{T}(a0))
        else
            nothing
        end
    end
    return GeometricGradientWorkspace{
        T, typeof(a_func), typeof(jac_fn), typeof(jac_cfg), typeof(linear_cache), typeof(a_const_lu),
    }(
        similar(path),
        zeros(T, Nx, N),
        zeros(T, N),
        zeros(T, N),
        zeros(T, Nx, N - 2),
        zeros(T, Nx, N - 2),
        zeros(T, Nx, N - 2),
        range(zero(T), one(T); length = N),
        a_func, a_const_lu, jac_fn, jac_cfg, J_buf, x_buf, linear_cache,
        zeros(T, Nx), zeros(T, Nx), zeros(T, Nx), zeros(T, Nx),
        zeros(T, Nx), zeros(T, Nx),
        zeros(T, Nx), zeros(T, Nx, Nx),
    )
end

# Allocation-free λ, θ computation. Buffers `Ainv_b`, `Ainv_φp` live in the workspace.
# `F_cached` may be a precomputed LU (constant non-diagonal a) or `nothing`.
function _lambda_theta!(θ_out, a_i, b_i, φp, Ainv_b, Ainv_φp, F_cached = nothing)
    if a_i isa LinearAlgebra.Diagonal
        d = a_i.diag
        @inbounds for k in eachindex(Ainv_b)
            Ainv_b[k] = b_i[k] / d[k]
            Ainv_φp[k] = φp[k] / d[k]
        end
    else
        F = F_cached === nothing ?
            LinearAlgebra.lu(a_i isa Matrix ? a_i : Matrix(a_i)) :
            F_cached
        copyto!(Ainv_b, b_i); LinearAlgebra.ldiv!(F, Ainv_b)
        copyto!(Ainv_φp, φp);  LinearAlgebra.ldiv!(F, Ainv_φp)
    end
    num = dot(b_i, Ainv_b)
    den = dot(φp, Ainv_φp)
    λ = den > 1.0e-28 ? sqrt(num / den) : zero(eltype(b_i))
    @inbounds for k in eachindex(θ_out)
        θ_out[k] = λ * Ainv_φp[k] - Ainv_b[k]
    end
    return λ
end

function _ws_jacobian!(ws, xi)
    copyto!(ws.x_buf, xi)
    ForwardDiff.jacobian!(ws.J_buf, ws.jac_fn, ws.x_buf, ws.jac_cfg, Val{false}())
    return ws.J_buf
end

function _assemble_explicit_rhs!(::Val{true}, ws, a_i, _θ, _b_i, φp, _xi, i)
    J = _ws_jacobian!(ws, _xi)
    LinearAlgebra.mul!(ws.Hpphi, J, φp)
    LinearAlgebra.mul!(ws.Hphi, J', _θ)
    LinearAlgebra.mul!(ws.aHphi, a_i, ws.Hphi)
    @inbounds for k in eachindex(ws.aHphi)
        ws.rhs_explicit[k, i - 1] = -ws.lambdas[i] * ws.Hpphi[k] + ws.aHphi[k]
    end
    return nothing
end

function _assemble_explicit_rhs!(::Val{false}, ws, a_i, θ, _b_i, φp, xi, i)
    J = _ws_jacobian!(ws, xi)
    LinearAlgebra.mul!(ws.Hpphi, J, φp)
    LinearAlgebra.mul!(ws.Hphi, J', θ)
    h = _fd_step(eltype(xi))
    inv_2h = 1 / (2 * h)
    Nx = length(ws.Hphi)
    xp = ws.x_probe; dla = ws.dla_buf
    copyto!(xp, xi)
    @inbounds for l in 1:Nx
        xp[l] = xi[l] + h
        a_xp = ws.a_func(xp)
        xp[l] = xi[l] - h
        a_xm = ws.a_func(xp)
        xp[l] = xi[l]
        @. dla = (a_xp - a_xm) * inv_2h
        ws.Hphi[l] += 0.5 * dot(θ, dla, θ)
        LinearAlgebra.mul!(ws.Hpphi, dla, θ, φp[l], one(eltype(ws.Hpphi)))
    end
    LinearAlgebra.mul!(ws.aHphi, a_i, ws.Hphi)
    @inbounds for k in eachindex(ws.aHphi)
        ws.rhs_explicit[k, i - 1] = -ws.lambdas[i] * ws.Hpphi[k] + ws.aHphi[k]
    end
    return nothing
end

"""
$(TYPEDSIGNATURES)

In-place one-step projected-gradient update of `ws.update` from `path`. Solves eq. (6) of
[heymann_pathways_2008](@citet) for an `N`-point path.

## Keyword arguments
  - `stepsize = 0.1`
  - `diff_order = 4`: order of the central finite-difference stencil along the path (2 or 4).
"""
function geometric_gradient_step!(
        ws::GeometricGradientWorkspace, sys, path; stepsize = 0.1, diff_order = 4,
    )
    N = size(path, 2); Nx = size(path, 1)
    dx = one(eltype(path)) / (N - 1)
    path_velocity!(ws.x_prime, path, ws.arc; order = diff_order)
    @views for i in 2:(N - 1)
        ws.drift_cache[:, i - 1] .= drift(sys, path[:, i])
    end
    @views for i in 2:(N - 1)
        xi = path[:, i]
        a_i = ws.a_func(xi)
        b_i = ws.drift_cache[:, i - 1]
        φp = ws.x_prime[:, i]
        ws.lambdas[i] = _lambda_theta!(
            view(ws.theta_cache, :, i - 1), a_i, b_i, φp, ws.Ainv_b, ws.Ainv_φp, ws.a_const_lu,
        )
    end
    @views for i in 2:(N - 1)
        ws.lambdas_prime[i] = (ws.lambdas[i + 1] - ws.lambdas[i - 1]) / (2 * dx)
    end
    @views for i in 2:(N - 1)
        xi = path[:, i]
        a_i = ws.a_func(xi)
        b_i = ws.drift_cache[:, i - 1]
        φp = ws.x_prime[:, i]
        θ = ws.theta_cache[:, i - 1]
        _assemble_explicit_rhs!(is_constant(ws.a_func), ws, a_i, θ, b_i, φp, xi, i)
        @inbounds for k in 1:Nx
            ws.rhs_explicit[k, i - 1] += ws.lambdas[i] * ws.lambdas_prime[i] * φp[k]
        end
    end
    _gmam_implicit_shared!(ws, path, N, stepsize, dx)
    return ws.update
end

function _gmam_implicit_shared!(ws::GeometricGradientWorkspace, path, N, stepsize, dx)
    Tmat = ws.linear_cache.A
    rhs = ws.linear_cache.b
    fill!(Tmat.d, 1); fill!(Tmat.dl, 0); fill!(Tmat.du, 0)
    @inbounds for i in 2:(N - 1)
        α = stepsize * ws.lambdas[i]^2 / dx^2
        if isfinite(α)
            Tmat.d[i] += 2α
            Tmat.dl[i - 1] = -α
            Tmat.du[i] = -α
        end
    end
    @inbounds for j in 1:size(path, 1)
        rhs[1] = path[j, 1]
        rhs[end] = path[j, end]
        for i in 2:(N - 1)
            rhs_val = path[j, i] + stepsize * ws.rhs_explicit[j, i - 1]
            rhs[i] = isfinite(rhs_val) ? rhs_val : path[j, i]
        end
        LinearSolve.reinit!(ws.linear_cache; A = Tmat, b = rhs)
        solve!(ws.linear_cache)
        @views ws.update[j, :] .= ws.linear_cache.u
    end
    return nothing
end
