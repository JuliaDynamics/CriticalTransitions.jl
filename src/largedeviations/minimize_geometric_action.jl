"""
    minimize_geometric_action(sys::CoupledSDEs, x_i, x_f, optimizer=GeometricGradient(); kwargs...)

Computes the minimizer of the geometric Freidlin-Wentzell action via the geometric minimum
action method (gMAM), using optimizers of Optimization.jl or the projected gradient method
of [heymann_pathways_2008](@citet) with optional backtracking.

The minimizer is computed for system `sys` over all paths from `x_i` to `x_f`. To set an
initial path different from a straight line, see the multiple-dispatch method
`minimize_geometric_action(sys::CoupledSDEs, init::Matrix, optimizer; kwargs...)`.

## Keyword arguments
  - `npoints::Int = 100`
  - `maxiters::Int = 100`
  - `abstol::Real = NaN`, `reltol::Real = NaN`
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

function minimize_geometric_action(sys::CoupledSDEs, init::Matrix; kwargs...)
    return minimize_geometric_action(sys, init, GeometricGradient(); kwargs...)
end

function minimize_geometric_action(
        sys::CoupledSDEs, init::Matrix, optimizer::GeometricGradient;
        maxiters::Int = 100, abstol::Real = NaN, reltol::Real = NaN,
        verbose = false, show_progress = true,
    )
    proper_FW_system(sys)
    path = deepcopy(init)
    N = size(init, 2)
    alpha = zeros(N)
    arc = range(0, 1.0; length = N)

    ws = geometric_gradient_workspace(sys, path)
    path_prev = similar(path)

    A_at = _action_metric(sys)
    b_fn = let sys = sys
        x -> drift(sys, x)
    end
    x_i = init[:, 1]; x_f = init[:, end]
    S(x) = _geometric_action_from_drift(b_fn, fix_ends(x, x_i, x_f), 1.0, A_at)

    interpolate_path!(path, alpha, arc)
    initial_action = S(path)

    function try_step!(ϵ)
        geometric_gradient_step!(ws, sys, path; stepsize = ϵ)
        path .= ws.update
        interpolate_path!(path, alpha, arc)
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
        sys::CoupledSDEs, init::Matrix, optimizer;
        maxiters::Int = 100, abstol::Real = NaN, reltol::Real = NaN,
        ad_type = OptimizationBase.AutoFiniteDiff(),
        show_progress = true,
    )
    proper_FW_system(sys)
    N = size(init, 2)
    alpha = zeros(N)
    arc = range(0, 1.0; length = N)
    A_at = _action_metric(sys)
    b_fn = let sys = sys
        x -> drift(sys, x)
    end
    x_i = init[:, 1]; x_f = init[:, end]
    S(x) = _geometric_action_from_drift(b_fn, fix_ends(x, x_i, x_f), 1.0, A_at)

    optf = SciMLBase.OptimizationFunction((x, _) -> S(x), ad_type)
    prob = SciMLBase.OptimizationProblem(optf, init, ())
    prog = Progress(maxiters; enabled = show_progress)
    function callback(state, _)
        interpolate_path!(state.u, alpha, arc)
        show_progress ? next!(prog) : nothing
        return false
    end
    sol = solve(prob, optimizer; maxiters, callback, abstol, reltol)
    path = sol.u
    return MinimumActionPath(StateSpaceSet(path'), S(path))
end

"""
$(TYPEDSIGNATURES)

Solves eq. (6) of [heymann_pathways_2008](@citet) for an initial `path` of `N` points.

## Keyword arguments
  - `stepsize = 0.1`
  - `diff_order = 4`
"""
function geometric_gradient_step(sys::CoupledSDEs, path, N; stepsize = 0.1, diff_order = 4)
    ws = geometric_gradient_workspace(sys, path)
    geometric_gradient_step!(ws, sys, path; stepsize = stepsize, diff_order = diff_order)
    return ws.update
end

struct GeometricGradientWorkspace{T, A, J, LC}
    update::Matrix{T}
    x_prime::Matrix{T}
    lambdas::Vector{T}
    lambdas_prime::Vector{T}
    drift_cache::Matrix{T}
    theta_cache::Matrix{T}
    rhs_explicit::Matrix{T}
    arc::StepRangeLen{T}
    a_func::A
    jac::J
    linear_cache::LC
    fd_probe::Vector{T}
    Hphi::Vector{T}
    Hpphi::Vector{T}
    aHphi::Vector{T}
end

function geometric_gradient_workspace(sys::CoupledSDEs, path)
    T = eltype(path)
    Nx, N = size(path)
    proper_FW_system(sys)
    a_func = _trace_normalized_a(sys)
    x_ref = collect(view(path, :, 1))
    _validate_and_classify_a(a_func, x_ref)
    params = current_parameters(sys)
    jac_fun = jacobian(sys)
    jac = let jac_fun = jac_fun, params = params
        x -> jac_fun(x, params, 0.0)
    end
    dl = zeros(T, N - 1); d = ones(T, N); du = zeros(T, N - 1)
    Tmat = LinearAlgebra.Tridiagonal(dl, d, du)
    rhs = zeros(T, N)
    linear_cache = init(
        LinearProblem(Tmat, rhs), LUFactorization();
        alias = SciMLBase.LinearAliasSpecifier(; alias_A = true, alias_b = true),
    )
    return GeometricGradientWorkspace{T, typeof(a_func), typeof(jac), typeof(linear_cache)}(
        similar(path),
        zeros(T, Nx, N),
        zeros(T, N),
        zeros(T, N),
        zeros(T, Nx, N - 2),
        zeros(T, Nx, N - 2),
        zeros(T, Nx, N - 2),
        range(zero(T), one(T); length = N),
        a_func, jac, linear_cache,
        zeros(T, Nx), zeros(T, Nx), zeros(T, Nx), zeros(T, Nx),
    )
end

function _lambda_theta!(θ_out, a_i, b_i, φp)
    F = LinearAlgebra.factorize(a_i isa LinearAlgebra.Diagonal ? a_i : Matrix(a_i))
    inv_b  = F \ b_i
    inv_φp = F \ φp
    num = dot(b_i, inv_b)
    den = dot(φp,  inv_φp)
    λ = den > 1.0e-28 ? sqrt(num / den) : zero(eltype(b_i))
    @inbounds for k in eachindex(θ_out)
        θ_out[k] = λ * inv_φp[k] - inv_b[k]
    end
    return λ
end

function _assemble_explicit_rhs!(::Val{true}, ws, a_i, θ, b_i, φp, xi, i)
    J = ws.jac(xi)
    LinearAlgebra.mul!(ws.Hpphi, J, φp)
    LinearAlgebra.mul!(ws.Hphi,  J', θ)
    LinearAlgebra.mul!(ws.aHphi, a_i, ws.Hphi)
    @inbounds for k in eachindex(ws.aHphi)
        ws.rhs_explicit[k, i - 1] = -ws.lambdas[i] * ws.Hpphi[k] + ws.aHphi[k]
    end
    return nothing
end

function _assemble_explicit_rhs!(::Val{false}, ws, a_i, θ, b_i, φp, xi, i)
    J = ws.jac(xi)
    LinearAlgebra.mul!(ws.Hpphi, J, φp)
    LinearAlgebra.mul!(ws.Hphi,  J', θ)
    h = _fd_step(eltype(xi))
    fill!(ws.fd_probe, 0)
    @inbounds for l in eachindex(ws.Hphi)
        ws.fd_probe[l] = h
        dla = (ws.a_func(xi .+ ws.fd_probe) .- ws.a_func(xi .- ws.fd_probe)) ./ (2 * h)
        ws.Hphi[l] += 0.5 * dot(θ, dla, θ)
        for k in eachindex(ws.Hpphi)
            for k2 in eachindex(θ)
                ws.Hpphi[k] += dla[k, k2] * θ[k2] * φp[l]
            end
        end
        ws.fd_probe[l] = 0
    end
    LinearAlgebra.mul!(ws.aHphi, a_i, ws.Hphi)
    @inbounds for k in eachindex(ws.aHphi)
        ws.rhs_explicit[k, i - 1] = -ws.lambdas[i] * ws.Hpphi[k] + ws.aHphi[k]
    end
    return nothing
end

function geometric_gradient_step!(
        ws::GeometricGradientWorkspace, sys, path; stepsize = 0.1, diff_order = 4,
    )
    N = size(path, 2)
    Nx = size(path, 1)
    dx = one(eltype(path)) / (N - 1)
    path_velocity!(ws.x_prime, path, ws.arc; order = diff_order)
    @views for i in 2:(N - 1)
        ws.drift_cache[:, i - 1] .= drift(sys, path[:, i])
    end
    @views for i in 2:(N - 1)
        xi  = path[:, i]
        a_i = ws.a_func(xi)
        b_i = ws.drift_cache[:, i - 1]
        φp  = ws.x_prime[:, i]
        ws.lambdas[i] = _lambda_theta!(view(ws.theta_cache, :, i - 1), a_i, b_i, φp)
    end
    @views for i in 2:(N - 1)
        ws.lambdas_prime[i] = (ws.lambdas[i + 1] - ws.lambdas[i - 1]) / (2 * dx)
    end
    @views for i in 2:(N - 1)
        xi  = path[:, i]
        a_i = ws.a_func(xi)
        b_i = ws.drift_cache[:, i - 1]
        φp  = ws.x_prime[:, i]
        θ   = ws.theta_cache[:, i - 1]
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
    rhs  = ws.linear_cache.b
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
        rhs[1]   = path[j, 1]
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
