"""
$(TYPEDSIGNATURES)

Computes the minimizer of the geometric Freidlin-Wentzell action based on the geometric
minimum action method (gMAM), using optimizers of OPtimization.jl or the original formulation
by Heymann and Vanden-Eijnden[^1]. Only Freidlin-Wentzell action has a geometric formulation.

To set an initial path different from a straight line, see the multiple dispatch method

  - `geometric_min_action_method(sys::CoupledSDEs, init::Matrix, arclength::Real; kwargs...)`.

## Keyword arguments

  - `maxiter::Int=100`: maximum number of optimization iterations before the algorithm stops
  - `abstol=1e-8`: absolute tolerance of action gradient to determine convergence
  - `reltol=1e-8`: relative tolerance of action gradient to determine convergence
  - `method = Adam()`: minimization algorithm (see below)
  - `=0.1`: step size parameter in gradient descent HeymannVandenEijnden implementation.
  - `verbose=false`: if true, print additional output
  - `show_progress=true`: if true, display a progress bar

## Minimization algorithms

The `method` keyword argument takes solver methods of the
[`Optimization.jl`](https://docs.sciml.ai/Optimization/) package; alternatively,
the option `method = "HeymannVandenEijnden"` implements the original gMAM
algorithm [heymann_pathways_2008](@cite).
"""
function geometric_min_action_method(
    sys::ContinuousTimeDynamicalSystem, x_i, x_f; N=100, kwargs...
)
    path = reduce(hcat, range(x_i, x_f; length=N))
    return geometric_min_action_method(sys, path; kwargs...)
end

"""
$(TYPEDSIGNATURES)

Runs the geometric Minimum Action Method (gMAM) to find the minimum action path (instanton) from an
initial condition `init`, given a system `sys` and total arc length `arclength`.

The initial path `init` must be a matrix of size `(D, N)`, where `D` is the dimension of the
system and `N` is the number of path points.

For more information see the main method,
[`geometric_min_action_method(sys::CoupledSDEs, x_i, x_f, arclength::Real; kwargs...)`](@ref).
"""
function geometric_min_action_method(
    sys::ContinuousTimeDynamicalSystem,
    init::Matrix;
    maxiter::Int=100,
    abstol=nothing,
    reltol=nothing,
    method=Optimisers.Adam(),
    AD=Optimization.AutoFiniteDiff(),
    ϵ=0.1,
    verbose=false,
    show_progress=true,
)
    if sys isa CoupledSDEs
        proper_MAM_system(sys)
    end

    path = deepcopy(init)
    N = length(init[1, :])
    alpha = zeros(N)
    arc = range(0, 1.0; length=N)

    S(x) = geometric_action(sys, fix_ends(x, init[:, 1], init[:, end]), 1.0)

    prog = Progress(maxiter; enabled=show_progress)
    if method == "HeymannVandenEijnden"
        for i in 1:maxiter
            update = heymann_vandeneijnden_step(sys, path, N; tau=ϵ)
            path .= update
            interpolate_path!(path, alpha, arc)
            next!(prog)
        end
    else
        optf = SciMLBase.OptimizationFunction((x, _) -> S(x), AD)
        prob = SciMLBase.OptimizationProblem(optf, init, ())

        function callback(state, loss_val)
            interpolate_path!(state.u, alpha, arc)
            show_progress ? next!(prog) : nothing
            return false
        end

        sol = solve(
            prob, method; maxiters=maxiter, callback=callback, abstol=abstol, reltol=reltol
        )
        path = sol.u
    end
    # if verbose && !converged
    #     @warn("Stopped after reaching maximum number of $(maxiter) iterations.")
    # end
    return MinimumActionPath(StateSpaceSet(path'), S(path))
end

"""
$(TYPEDSIGNATURES)

Solves eq. (6) of Ref.[^1] for an initial `path` with `N` points and arclength `L`.

## Keyword arguments

  - `tau = 0.1`: step size
  - `diff_order = 4`: order of the finite differencing along the path. Either `2` or `4`.

[^1]: [Heymann and Vanden-Eijnden, PRL (2008)](https://link.aps.org/doi/10.1103/PhysRevLett.100.140601)
"""
function heymann_vandeneijnden_step(
    sys::ContinuousTimeDynamicalSystem, path, N; tau=0.1, diff_order=4
)
    L = 1.0
    dx = L / (N - 1)
    update = zeros(size(path))
    lambdas, lambdas_prime = zeros(N), zeros(N)
    x_prime = path_velocity(path, range(0, L; length=N); order=diff_order)

    A = inv(covariance_matrix(sys))

    for i in 2:(N - 1)
        velocity_norm = anorm(x_prime[:, i], A)
        if velocity_norm > 1e-14
            lambdas[i] = anorm(drift(sys, path[:, i]), A) / velocity_norm
            if !isfinite(lambdas[i])
                lambdas[i] = 0.0
            end
        end
    end
    for i in 2:(N - 1)
        lambdas_prime[i] = (lambdas[i + 1] - lambdas[i - 1]) / (2 * dx)
    end

    b(x) = drift(sys, x)
    params = current_parameters(sys)
    jac_fun = jacobian(sys)

    jac = let jac_fun=jac_fun, params=params
        x -> jac_fun(x, params, 0.0)
    end

    J = [jac(path[:, i]) for i in 2:(N - 1)]
    prod1 = [(J[i - 1] - J[i - 1]') * x_prime[:, i] for i in 2:(N - 1)]
    prod2 = [(J[i - 1]') * b(path[:, i]) for i in 2:(N - 1)]

    # Solve linear system M*x = v for each system dimension
    #! might be made faster using LinearSolve.jl special solvers
    for j in 1:size(path, 1)
        M = Matrix{Float64}(I(N))
        v = zeros(N)

        # Boundary conditions
        v[1] = path[j, 1]
        v[end] = path[j, end]

        # Linear system of equations
        for i in 2:(N - 1)
            alpha = tau * lambdas[i]^2 / (dx^2)
            if !isfinite(alpha)
                v[i] = path[j, i]
                continue
            end
            M[i, i] += 2 * alpha
            M[i, i - 1] = -alpha
            M[i, i + 1] = -alpha

            rhs =
                path[j, i] +
                tau * (
                    -lambdas[i] * prod1[i - 1][j] - prod2[i - 1][j] +
                    lambdas[i] * lambdas_prime[i] * x_prime[j, i]
                )
            v[i] = isfinite(rhs) ? rhs : path[j, i]
        end

        # Solve and store solution
        sol = M \ v
        update[j, :] = sol
    end
    return update
end
