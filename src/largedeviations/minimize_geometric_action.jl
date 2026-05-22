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
        proper_FW_system(sys)
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

struct GeometricGradientWorkspace{
        Tupdate, Tprime, Tvec, Tmat, Ttmp, Tarc, TA, Tjac, AF, NS <: NoiseShape,
    }
    update::Tupdate
    x_prime::Tprime
    lambdas::Tvec
    lambdas_prime::Tvec
    prod1::Tmat
    prod2::Tmat
    drift_cache::Tmat
    theta_cache::Tmat
    rhs_explicit::Tmat
    temp1::Ttmp
    temp2::Ttmp
    alpha_cache::Tvec
    dl::Tvec
    d::Tvec
    du::Tvec
    v::Tvec
    arc::Tarc
    A::TA
    jac::Tjac
    a_func::AF
end

# Build a trace-normalized `a(x)` callable for the gMAM workspace. Mirrors the logic
# in FreidlinWentzellHamiltonian's constructor: identity for CoupledODEs, constant
# wrapped in `Returns` for additive SDE, closure for state-dependent SDE.
function _build_a_func(sys::ContinuousTimeDynamicalSystem)
    D = dimension(sys)
    if sys isa CoupledODEs
        return Returns(LinearAlgebra.Diagonal(ones(Float64, D)))
    end
    σ_fn = diffusion_function(sys)
    ps   = current_parameters(sys)
    u₀   = current_state(sys)
    σ0 = σ_fn(u₀, ps, 0.0)
    σ_mat0 = σ0 isa AbstractMatrix ? σ0 : LinearAlgebra.Diagonal(σ0)
    a0 = σ_mat0 * σ_mat0'
    s = LinearAlgebra.tr(a0) / D
    if sys.noise_type[:additive]
        a_const = LinearAlgebra.isdiag(a0) ?
            LinearAlgebra.Diagonal(collect(LinearAlgebra.diag(a0)) ./ s) :
            Matrix(a0 ./ s)
        return Returns(a_const)
    else
        return let σ_fn = σ_fn, ps = ps, s = s
            x -> begin
                σx = σ_fn(x, ps, 0.0)
                σ_mat = σx isa AbstractMatrix ? σx : LinearAlgebra.Diagonal(σx)
                (σ_mat * σ_mat') / s
            end
        end
    end
end

function geometric_gradient_workspace(sys::ContinuousTimeDynamicalSystem, path)
    N = size(path, 2)
    NS = _classify_noise_shape(sys)
    a_func = _build_a_func(sys)
    A = a_func(view(path, :, 1))
    A = LinearAlgebra.isdiag(A) && size(A, 1) == size(A, 2) ?
        inv(LinearAlgebra.Diagonal(LinearAlgebra.diag(A))) : inv(Matrix(A))
    params = current_parameters(sys)
    jac_fun = jacobian(sys)
    jac = let jac_fun = jac_fun, params = params
        x -> jac_fun(x, params, 0.0)
    end

    update = similar(path)
    x_prime = zeros(size(path))
    lambdas = zeros(N)
    lambdas_prime = zeros(N)
    prod1 = zeros(size(path, 1), N - 2)
    prod2 = zeros(size(path, 1), N - 2)
    drift_cache = zeros(size(path, 1), N - 2)
    theta_cache = zeros(size(path, 1), N - 2)
    rhs_explicit = zeros(size(path, 1), N - 2)
    temp1 = zeros(size(path, 1))
    temp2 = zeros(size(path, 1))
    alpha_cache = zeros(N)
    dl = zeros(N - 1)
    d = ones(N)
    du = zeros(N - 1)
    v = zeros(N)
    arc = range(0, 1.0; length = N)

    return GeometricGradientWorkspace{
        typeof(update), typeof(x_prime), typeof(lambdas), typeof(prod1),
        typeof(temp1), typeof(arc), typeof(A), typeof(jac), typeof(a_func), typeof(NS),
    }(
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
        ws::GeometricGradientWorkspace{Tupdate, Tprime, Tvec, Tmat, Ttmp, Tarc, TA, Tjac, AF, NS},
        sys, path;
        stepsize = 0.1, diff_order = 4,
    ) where {Tupdate, Tprime, Tvec, Tmat, Ttmp, Tarc, TA, Tjac, AF, NS <: NoiseShape}
    _geometric_gradient_step!(ws, sys, path, NS(); stepsize, diff_order)
    return ws.update
end

# Additive branch: original main-branch logic, using the cached constant A = inv(a).
function _geometric_gradient_step!(
        ws, sys, path, ::AdditiveNoise; stepsize = 0.1, diff_order = 4,
    )
    N = size(path, 2)
    dx = 1.0 / (N - 1)

    path_velocity!(ws.x_prime, path, ws.arc; order = diff_order)

    @views for i in 2:(N - 1)
        velocity_norm = anorm(ws.x_prime[:, i], ws.A)
        if velocity_norm > 1.0e-14
            ws.drift_cache[:, i - 1] .= drift(sys, path[:, i])
            ws.lambdas[i] = anorm(ws.drift_cache[:, i - 1], ws.A) / velocity_norm
            isfinite(ws.lambdas[i]) || (ws.lambdas[i] = 0.0)
        else
            ws.lambdas[i] = 0.0
        end
    end
    @views for i in 2:(N - 1)
        ws.lambdas_prime[i] = (ws.lambdas[i + 1] - ws.lambdas[i - 1]) / (2 * dx)
    end
    @views for i in 2:(N - 1)
        J = ws.jac(path[:, i])
        LinearAlgebra.mul!(ws.temp1, J, ws.x_prime[:, i])
        LinearAlgebra.mul!(ws.temp2, J', ws.x_prime[:, i])
        ws.prod1[:, i - 1] .= ws.temp1 .- ws.temp2
        LinearAlgebra.mul!(ws.prod2[:, i - 1], J', ws.drift_cache[:, i - 1])
    end

    _gmam_implicit_shared!(ws, path, N, stepsize, dx; is_state_dep = false)
    return nothing
end

# DiagonalNoise: per-grid-point diagonal a(x_i). θ_i recovered via diag(a(x_i)),
# explicit RHS includes (1/2) ∂_l a_kk θ_k² and (∂_l a_kk θ_k) φ'_l contributions.
function _geometric_gradient_step!(
        ws, sys, path, ::DiagonalNoise; stepsize = 0.1, diff_order = 4,
    )
    N  = size(path, 2)
    Nx = size(path, 1)
    dx = 1.0 / (N - 1)
    h_fd = max(sqrt(eps(eltype(path))), eltype(path)(1.0e-8))

    path_velocity!(ws.x_prime, path, ws.arc; order = diff_order)

    @views for i in 2:(N - 1)
        ws.drift_cache[:, i - 1] .= drift(sys, path[:, i])
    end

    @views for i in 2:(N - 1)
        xi = path[:, i]
        a_i = LinearAlgebra.diag(ws.a_func(xi))
        b_i = ws.drift_cache[:, i - 1]
        φp  = ws.x_prime[:, i]
        num = sum(b_i .^ 2 ./ a_i)
        den = sum(φp .^ 2 ./ a_i)
        λ = den > 1.0e-28 ? sqrt(num / den) : 0.0
        isfinite(λ) || (λ = 0.0)
        ws.lambdas[i] = λ
        @inbounds for k in 1:Nx
            ws.theta_cache[k, i - 1] = (λ * φp[k] - b_i[k]) / a_i[k]
        end
    end
    @views for i in 2:(N - 1)
        ws.lambdas_prime[i] = (ws.lambdas[i + 1] - ws.lambdas[i - 1]) / (2 * dx)
    end

    e = zeros(eltype(path), Nx)
    Hphi  = zeros(eltype(path), Nx)
    Hpphi = zeros(eltype(path), Nx)
    @views for i in 2:(N - 1)
        xi = path[:, i]
        φp = ws.x_prime[:, i]
        θ  = ws.theta_cache[:, i - 1]
        J  = ws.jac(xi)
        fill!(Hphi, 0); fill!(Hpphi, 0)
        LinearAlgebra.mul!(Hpphi, J, φp)
        LinearAlgebra.mul!(Hphi,  J', θ)
        @inbounds for l in 1:Nx
            fill!(e, 0); e[l] = h_fd
            ap = LinearAlgebra.diag(ws.a_func(xi .+ e))
            am = LinearAlgebra.diag(ws.a_func(xi .- e))
            @inbounds for k in 1:Nx
                dla = (ap[k] - am[k]) / (2 * h_fd)
                Hphi[l]  += 0.5 * dla * θ[k]^2
                Hpphi[k] += dla * θ[k] * φp[l]
            end
        end
        a_i = LinearAlgebra.diag(ws.a_func(xi))
        λ  = ws.lambdas[i]
        λp = ws.lambdas_prime[i]
        @inbounds for k in 1:Nx
            ws.rhs_explicit[k, i - 1] = -λ * Hpphi[k] + a_i[k] * Hphi[k] + λ * λp * φp[k]
        end
    end

    _gmam_implicit_shared!(ws, path, N, stepsize, dx; is_state_dep = true)
    return nothing
end

# GeneralNoise: per-grid-point matrix a(x_i). θ_i = a(x_i)^{-1} (λ φ' - b).
# Explicit RHS uses (1/2) θ^T (∂_l a) θ for each l, plus (∂_l a · θ) φ'_l terms.
function _geometric_gradient_step!(
        ws, sys, path, ::GeneralNoise; stepsize = 0.1, diff_order = 4,
    )
    N  = size(path, 2)
    Nx = size(path, 1)
    dx = 1.0 / (N - 1)
    h_fd = max(sqrt(eps(eltype(path))), eltype(path)(1.0e-8))

    path_velocity!(ws.x_prime, path, ws.arc; order = diff_order)

    @views for i in 2:(N - 1)
        ws.drift_cache[:, i - 1] .= drift(sys, path[:, i])
    end

    @views for i in 2:(N - 1)
        xi = path[:, i]
        a_i = ws.a_func(xi)
        b_i = ws.drift_cache[:, i - 1]
        φp  = ws.x_prime[:, i]
        A_inv_b   = a_i \ Vector(b_i)
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
    @views for i in 2:(N - 1)
        ws.lambdas_prime[i] = (ws.lambdas[i + 1] - ws.lambdas[i - 1]) / (2 * dx)
    end

    e = zeros(eltype(path), Nx)
    Hphi  = zeros(eltype(path), Nx)
    Hpphi = zeros(eltype(path), Nx)
    aHphi = zeros(eltype(path), Nx)
    @views for i in 2:(N - 1)
        xi = path[:, i]
        φp = ws.x_prime[:, i]
        θ  = ws.theta_cache[:, i - 1]
        J  = ws.jac(xi)
        a_i = ws.a_func(xi)
        fill!(Hphi, 0); fill!(Hpphi, 0)
        LinearAlgebra.mul!(Hpphi, J, φp)
        LinearAlgebra.mul!(Hphi,  J', θ)
        @inbounds for l in 1:Nx
            fill!(e, 0); e[l] = h_fd
            ap = ws.a_func(xi .+ e)
            am = ws.a_func(xi .- e)
            @inbounds for k1 in 1:Nx, k2 in 1:Nx
                dla = (ap[k1, k2] - am[k1, k2]) / (2 * h_fd)
                Hphi[l]   += 0.5 * dla * θ[k1] * θ[k2]
                Hpphi[k1] += dla * θ[k2] * φp[l]
            end
        end
        LinearAlgebra.mul!(aHphi, a_i, Hphi)
        λ  = ws.lambdas[i]; λp = ws.lambdas_prime[i]
        @inbounds for k in 1:Nx
            ws.rhs_explicit[k, i - 1] = -λ * Hpphi[k] + aHphi[k] + λ * λp * φp[k]
        end
    end

    _gmam_implicit_shared!(ws, path, N, stepsize, dx; is_state_dep = true)
    return nothing
end

# Shared implicit step `(I - ε λ_i² ∂_s² / dx²) x_new = R` acting per DOF.
# Row i has α_i = ε λ_i²/dx², with both off-diagonals from row i using α_i (matching
# the main-branch convention). The `α_i == 0` skip preserves λ=0 path points.
function _gmam_implicit_shared!(ws, path, N, stepsize, dx; is_state_dep)
    ws.d  .= 1
    ws.dl .= 0
    ws.du .= 0
    @views for i in 2:(N - 1)
        ws.alpha_cache[i] = stepsize * ws.lambdas[i]^2 / dx^2
        if isfinite(ws.alpha_cache[i])
            ws.d[i]      += 2 * ws.alpha_cache[i]
            ws.dl[i - 1]  = -ws.alpha_cache[i]
            ws.du[i]      = -ws.alpha_cache[i]
        else
            ws.alpha_cache[i] = 0.0
        end
    end

    T = LinearAlgebra.Tridiagonal(ws.dl, ws.d, ws.du)
    prob = LinearProblem(T, ws.v)
    cache = init(
        prob, LUFactorization();
        alias = SciMLBase.LinearAliasSpecifier(; alias_A = true, alias_b = true),
    )
    rhs = cache.b

    @inbounds for j in 1:size(path, 1)
        rhs[1]   = path[j, 1]
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
