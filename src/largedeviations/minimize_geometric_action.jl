"""
    minimize_geometric_action(sys::CoupledSDEs, x_i, x_f, optimizer=GeometricGradient(); kwargs...)

Computes the minimizer of the geometric Freidlin-Wentzell action based on the geometric
minimum action method (gMAM), using optimizers of Optimization.jl or the original formulation
by Heymann and Vanden-Eijnden [heymann_pathways_2008](@citet), which performs a projected gradient
descent.

The minimizer is computed for the system `sys` over all paths leading from the initial state
`x_i` to the final state `x_f`.
To set an initial path different from a straight line, see the multiple dispatch method

  - `minimize_geometric_action(sys::CoupledSDEs, init::Matrix, arclength::Real; kwargs...)`.

Returns a [`MinimumActionPath`](@ref) object.

## Keyword arguments

  - `npoints::Int=100`: number of discretization points along the path
  - `maxiters::Int=100`: maximum number of optimization iterations before the algorithm stops
  - `abstol::Real=NaN`: absolute tolerance of action change to determine convergence
  - `reltol::Real=NaN`: relative tolerance of action change to determine convergence
  - `ad_type=OptimizationBase.AutoFiniteDiff()`: type of automatic differentiation (only used with Optimization.jl solvers)
  - `verbose=false`: if true, print additional output
  - `show_progress=true`: if true, display a progress bar

The optional positional argument `optimizer` selects the minimization algorithm. It defaults to
`GeometricGradient()`, which enables backtracking step-size control (see [`GeometricGradient`](@ref)).
Pass `GeometricGradient(; max_backtracks=0)` to use a fixed step size, or any
Optimization.jl optimizer to use Optimization.jl (`ad_type` is only used in that case).
When backtracking is enabled, prefer a **large** initial step size via
`GeometricGradient(; stepsize=...)`: rejected steps are cheap and the controller reduces
the step size automatically, so starting large gives fast early progress without
sacrificing accuracy.

## Minimization algorithms

The `optimizer` argument accepts:
  - `GeometricGradient()`: projected-gradient descent with backtracking [heymann_pathways_2008](@citet)
  - Any solver from [`Optimization.jl`](https://docs.sciml.ai/Optimization/) (e.g., `OptimizationOptimisers.Adam()`)

## Notes
Only the Freidlin-Wentzell action has a geometric formulation.
"""
function minimize_geometric_action(
        sys::ContinuousTimeDynamicalSystem,
        x_i,
        x_f,
        optimizer = GeometricGradient();
        npoints::Int = 100,
        kwargs...,
    )
    path = reduce(hcat, range(x_i, x_f; length = npoints))
    return minimize_geometric_action(sys, path, optimizer; kwargs...)
end

function _gmam_setup(sys, init, maxiters, show_progress)
    if sys isa CoupledSDEs
        proper_MAM_system(sys)
    end
    path = deepcopy(init)
    N = length(init[1, :])
    alpha = zeros(N)
    arc = range(0, 1.0; length = N)
    S(x) = geometric_action(sys, fix_ends(x, init[:, 1], init[:, end]), 1.0)
    prog = Progress(maxiters; enabled = show_progress)
    return path, N, alpha, arc, S, prog
end

function minimize_geometric_action(
        sys::ContinuousTimeDynamicalSystem, init::Matrix; kwargs...
    )
    return minimize_geometric_action(sys, init, GeometricGradient(); kwargs...)
end

"""
$(TYPEDSIGNATURES)

Runs the geometric Minimum Action Method (gMAM) to find the minimum action path (instanton) from an
initial condition `init`, given a system `sys` and total arc length `arclength`.

The initial path `init` must be a matrix of size `(D, N)`, where `D` is the dimension of the
system and `N` is the number of path points.

For more information see the main method,
[`minimize_geometric_action(sys::CoupledSDEs, x_i, x_f, arclength::Real; kwargs...)`](@ref).
"""
function minimize_geometric_action(
        sys::ContinuousTimeDynamicalSystem,
        init::Matrix,
        optimizer::GeometricGradient;
        maxiters::Int = 100,
        abstol::Real = NaN,
        reltol::Real = NaN,
        verbose = false,
        show_progress = true,
    )
    path, _, alpha, arc, S, _ = _gmam_setup(sys, init, maxiters, false)

    ws = geometric_gradient_workspace(sys, path)
    path_prev = similar(path)

    # Ensure a consistent starting path for action comparisons
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
        optimizer,
        try_step!,
        save!,
        restore!,
        initial_action;
        maxiters,
        abstol,
        reltol,
        verbose,
        show_progress,
    )
    return MinimumActionPath(StateSpaceSet(path'), current_action)
end

function minimize_geometric_action(
        sys::ContinuousTimeDynamicalSystem,
        init::Matrix,
        optimizer;
        maxiters::Int = 100,
        abstol::Real = NaN,
        reltol::Real = NaN,
        ad_type = OptimizationBase.AutoFiniteDiff(),
        show_progress = true,
    )
    path, N, alpha, arc, S, prog = _gmam_setup(sys, init, maxiters, show_progress)

    optf = SciMLBase.OptimizationFunction((x, _) -> S(x), ad_type)
    prob = SciMLBase.OptimizationProblem(optf, init, ())

    function callback(state, loss_val)
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

Solves eq. (6) of [heymann_pathways_2008](@citet) for an initial `path` with `N` points and
arclength `L`.

## Keyword arguments

  - `stepsize = 0.1`: step size for gradient descent
  - `diff_order = 4`: order of the finite differencing along the path. Either `2` or `4`.

## References

- [heymann_pathways_2008](@cite)
"""
function geometric_gradient_step(
        sys::ContinuousTimeDynamicalSystem, path, N; stepsize = 0.1, diff_order = 4
    )
    ws = geometric_gradient_workspace(sys, path)
    geometric_gradient_step!(ws, sys, path; stepsize = stepsize, diff_order = diff_order)
    return ws.update
end

struct GeometricGradientWorkspace{Tupdate, Tprime, Tvec, Tmat, Ttmp, Tarc, TA, Tjac, TAF}
    update::Tupdate
    x_prime::Tprime
    lambdas::Tvec
    lambdas_prime::Tvec
    prod1::Tmat
    prod2::Tmat
    drift_cache::Tmat
    theta_cache::Tmat        # state-dependent: ϑ_i = a^{-1}(x_i)(λ_i φ'_i - b_i)
    rhs_explicit::Tmat       # state-dependent: per-grid-point explicit RHS contribution
    temp1::Ttmp
    temp2::Ttmp
    alpha_cache::Tvec
    dl::Tvec
    d::Tvec
    du::Tvec
    v::Tvec
    arc::Tarc
    A::TA                    # additive: inverse normalized covariance; multiplicative: nothing
    jac::Tjac
    a_func::TAF              # nothing (additive), or x -> a(x) (Vector for diagonal, Matrix for general)
end

function geometric_gradient_workspace(sys::ContinuousTimeDynamicalSystem, path)
    N = size(path, 2)
    Nx = size(path, 1)
    params = current_parameters(sys)
    jac_fun = jacobian(sys)
    jac = let jac_fun = jac_fun, params = params
        x -> jac_fun(x, params, 0.0)
    end

    # Convention B (Freidlin-Wentzell rate function): both branches divide the diffusion
    # tensor by the noise-magnitude scale `s = L_1(a(x_0))/D` extracted at the reference
    # state. See `_normalized_a_shape_func` in utils.jl.
    if sys isa CoupledSDEs && !sys.noise_type[:additive]
        a_func, _ = _normalized_a_shape_func(sys)
        A = nothing
    else
        a = covariance_matrix(sys)
        s = LinearAlgebra.norm(a, 1) / size(a, 1)
        A = s * inv(a)
        a_func = nothing
    end

    update = similar(path)
    x_prime = zeros(size(path))
    lambdas = zeros(N)
    lambdas_prime = zeros(N)
    prod1 = zeros(Nx, N - 2)
    prod2 = zeros(Nx, N - 2)
    drift_cache = zeros(Nx, N - 2)
    theta_cache = zeros(Nx, N - 2)
    rhs_explicit = zeros(Nx, N - 2)
    temp1 = zeros(Nx)
    temp2 = zeros(Nx)
    alpha_cache = zeros(N)
    dl = zeros(N - 1)
    d = ones(N)
    du = zeros(N - 1)
    v = zeros(N)
    arc = range(0, 1.0; length = N)

    return GeometricGradientWorkspace(
        update,
        x_prime,
        lambdas,
        lambdas_prime,
        prod1,
        prod2,
        drift_cache,
        theta_cache,
        rhs_explicit,
        temp1,
        temp2,
        alpha_cache,
        dl,
        d,
        du,
        v,
        arc,
        A,
        jac,
        a_func,
    )
end

function geometric_gradient_step!(
        ws::GeometricGradientWorkspace,
        sys::ContinuousTimeDynamicalSystem,
        path;
        stepsize = 0.1,
        diff_order = 4,
    )
    # Heymann-VandenEijnden gMAM relaxation step. Eq. (19) of Grafke-Schäfer-Vanden-Eijnden
    # 2017 gives the descent direction
    #     ∂φ/∂τ = λ² φ'' - λ H_ϑφ φ' + H_ϑϑ H_φ + λ λ' φ'.
    # For additive a = I, this collapses to -λ(J-Jᵀ)φ' - Jᵀb + λλ'φ' on the RHS and
    # T = I - ε λ² ∂²_s on the LHS. For state-dependent a(x) we keep the same descent
    # equation but evaluate H_ϑφ and a·H_φ pointwise with finite-difference ∂_x a.
    N = size(path, 2)
    Nx = size(path, 1)
    dx = 1.0 / (N - 1)

    path_velocity!(ws.x_prime, path, ws.arc; order = diff_order)

    is_state_dep = ws.a_func !== nothing

    # Compute λ_i and (for state-dependent) ϑ_i at each interior grid point.
    @views for i in 2:(N - 1)
        ws.drift_cache[:, i - 1] .= drift(sys, path[:, i])
    end

    if is_state_dep
        _gmam_state_dep_lambda_theta!(ws, path, N, Nx)
    else
        @views for i in 2:(N - 1)
            velocity_norm = anorm(ws.x_prime[:, i], ws.A)
            if velocity_norm > 1.0e-14
                ws.lambdas[i] = anorm(ws.drift_cache[:, i - 1], ws.A) / velocity_norm
                isfinite(ws.lambdas[i]) || (ws.lambdas[i] = 0.0)
            else
                ws.lambdas[i] = 0.0
            end
        end
    end

    @views for i in 2:(N - 1)
        ws.lambdas_prime[i] = (ws.lambdas[i + 1] - ws.lambdas[i - 1]) / (2 * dx)
    end

    # Build explicit RHS contribution (everything except the implicit λ² φ'' term).
    if is_state_dep
        _gmam_state_dep_rhs!(ws, path, N, Nx)
    else
        @views for i in 2:(N - 1)
            J = ws.jac(path[:, i])
            LinearAlgebra.mul!(ws.temp1, J, ws.x_prime[:, i])
            LinearAlgebra.mul!(ws.temp2, J', ws.x_prime[:, i])
            ws.prod1[:, i - 1] .= ws.temp1 .- ws.temp2
            LinearAlgebra.mul!(ws.prod2[:, i - 1], J', ws.drift_cache[:, i - 1])
        end
    end

    # Implicit step. Per Heymann-Vanden-Eijnden 2008 eq. (3.6) and Grafke et al. 2017
    # eq. (19), the gMAM descent contains the term λ² φ'' (no a^{-1} factor), so the
    # implicit operator (1 - h λ² ∂_s²) is **shared across DOFs** regardless of whether
    # a(x) is state-dependent. State-dependence enters only through the explicit RHS
    # built above. (Contrast sgMAM eq. (43), where the implicit operator is
    # (1 - h μ^{-2} H_ϑϑ^{-1} ∂_s²) and does carry a^{-1}.)
    _gmam_implicit_shared!(ws, path, N, stepsize, dx, is_state_dep)

    return ws.update
end

# ----- state-dependent: λ_i and ϑ_i = a^{-1}(λ φ' - b) per interior grid point -----
function _gmam_state_dep_lambda_theta!(ws, path, N, Nx)
    for i in 2:(N - 1)
        a_i = ws.a_func(view(path, :, i))
        b_i = view(ws.drift_cache, :, i - 1)
        φp = view(ws.x_prime, :, i)
        if _is_diagonal_a(a_i)
            a_vec = _diag_view(a_i)
            num = sum(b_i .^ 2 ./ a_vec)
            den = sum(φp .^ 2 ./ a_vec)
            λ = den > 1.0e-28 ? sqrt(num / den) : 0.0
            isfinite(λ) || (λ = 0.0)
            ws.lambdas[i] = λ
            @inbounds for k in 1:Nx
                ws.theta_cache[k, i - 1] = (λ * φp[k] - b_i[k]) / a_vec[k]
            end
        else
            A_inv_b = a_i \ Vector(b_i)
            A_inv_phi = a_i \ Vector(φp)
            num = dot(b_i, A_inv_b)
            den = dot(φp, A_inv_phi)
            λ = den > 1.0e-28 ? sqrt(num / den) : 0.0
            isfinite(λ) || (λ = 0.0)
            ws.lambdas[i] = λ
            θ = a_i \ (λ .* Vector(φp) .- Vector(b_i))
            @inbounds for k in 1:Nx
                ws.theta_cache[k, i - 1] = θ[k]
            end
        end
    end
    return nothing
end

# ----- state-dependent: explicit RHS contribution per interior grid point -----
# For each i, store ws.rhs_explicit[:, i-1] = -λ_i H_ϑφ φ'_i + a(x_i) H_φ_i + λ_i λ'_i φ'_i.
function _gmam_state_dep_rhs!(ws, path, N, Nx)
    h = max(sqrt(eps(eltype(path))), 1.0e-8)
    e = zeros(eltype(path), Nx)
    Hphi = zeros(eltype(path), Nx)            # H_φ at current grid point
    Hpphi = zeros(eltype(path), Nx)           # H_ϑφ φ' at current grid point
    aHphi = zeros(eltype(path), Nx)           # a(x) · H_φ at current grid point

    for i in 2:(N - 1)
        xi = path[:, i]
        φp = view(ws.x_prime, :, i)
        ϑ = view(ws.theta_cache, :, i - 1)
        J = ws.jac(xi)

        a_i = ws.a_func(xi)
        is_diag = _is_diagonal_a(a_i)

        fill!(Hphi, 0)
        fill!(Hpphi, 0)
        # (J φ') and (J^T ϑ) contributions (always present)
        LinearAlgebra.mul!(Hpphi, J, φp)                  # (J φ')_j
        LinearAlgebra.mul!(Hphi, J', ϑ)                   # (J^T ϑ)_l

        # ∂_x a contributions via central finite differences. For each direction l in 1..Nx,
        # compute ∂_l a; accumulate (1/2) ⟨ϑ, ∂_l a ϑ⟩ into H_φ_l and (∂_l a · ϑ) φ'_l into
        # (H_ϑφ φ').
        @inbounds for l in 1:Nx
            fill!(e, 0)
            e[l] = h
            ap = ws.a_func(xi .+ e)
            am = ws.a_func(xi .- e)
            if is_diag
                # ap, am are diagonal: ∂_l a_k = (ap[k] - am[k])/(2h)
                ap_v = _diag_view(ap)
                am_v = _diag_view(am)
                @inbounds for k in 1:Nx
                    dla = (ap_v[k] - am_v[k]) / (2 * h)
                    Hphi[l] += 0.5 * dla * ϑ[k]^2
                    Hpphi[k] += dla * ϑ[k] * φp[l]
                end
            else
                # ap, am are full Matrices; ∂_l a is a D×D matrix
                @inbounds for k1 in 1:Nx, k2 in 1:Nx
                    dla = (ap[k1, k2] - am[k1, k2]) / (2 * h)
                    Hphi[l] += 0.5 * dla * ϑ[k1] * ϑ[k2]
                    Hpphi[k1] += dla * ϑ[k2] * φp[l]
                end
            end
        end

        # a(x_i) · H_φ
        if is_diag
            a_vec = _diag_view(a_i)
            @inbounds for k in 1:Nx
                aHphi[k] = a_vec[k] * Hphi[k]
            end
        else
            LinearAlgebra.mul!(aHphi, a_i, Hphi)
        end

        # RHS at point i: -λ H_ϑφ φ' + a H_φ + λ λ' φ'
        λ = ws.lambdas[i]
        λp = ws.lambdas_prime[i]
        @inbounds for k in 1:Nx
            ws.rhs_explicit[k, i - 1] = -λ * Hpphi[k] + aHphi[k] + λ * λp * φp[k]
        end
    end
    return nothing
end

# ----- implicit step: (1 - h λ² ∂_s²) — shared tridiagonal across DOFs (eq. 19 gMAM) -----
# Whether a(x) is additive or state-dependent, the gMAM implicit operator is independent
# of a; state dependence enters only through the explicit RHS (precomputed above). The
# additive branch keeps the fast (J - Jᵀ)φ' / Jᵀb decomposition; the state-dependent
# branch uses the precomputed `rhs_explicit` matrix.
function _gmam_implicit_shared!(ws, path, N, stepsize, dx, is_state_dep)
    ws.d .= 1
    ws.dl .= 0
    ws.du .= 0
    @views for i in 2:(N - 1)
        ws.alpha_cache[i] = stepsize * ws.lambdas[i]^2 / (dx^2)
        if isfinite(ws.alpha_cache[i])
            ws.d[i] += 2 * ws.alpha_cache[i]
            ws.dl[i - 1] = -ws.alpha_cache[i]
            ws.du[i] = -ws.alpha_cache[i]
        else
            ws.alpha_cache[i] = 0.0
        end
    end

    T = LinearAlgebra.Tridiagonal(ws.dl, ws.d, ws.du)
    prob = LinearProblem(T, ws.v)
    cache = init(
        prob,
        LUFactorization();
        alias = SciMLBase.LinearAliasSpecifier(; alias_A = true, alias_b = true),
    )
    rhs = cache.b

    @inbounds for j in 1:size(path, 1)
        rhs[1] = path[j, 1]
        rhs[end] = path[j, end]
        for i in 2:(N - 1)
            if ws.alpha_cache[i] == 0.0
                rhs[i] = path[j, i]
                continue
            end
            rhs_val = if is_state_dep
                path[j, i] + stepsize * ws.rhs_explicit[j, i - 1]
            else
                path[j, i] + stepsize * (
                    -ws.lambdas[i] * ws.prod1[j, i - 1] - ws.prod2[j, i - 1] +
                        ws.lambdas[i] * ws.lambdas_prime[i] * ws.x_prime[j, i]
                )
            end
            rhs[i] = isfinite(rhs_val) ? rhs_val : path[j, i]
        end
        solve!(cache)
        ws.update[j, :] .= cache.u
    end
    return nothing
end
