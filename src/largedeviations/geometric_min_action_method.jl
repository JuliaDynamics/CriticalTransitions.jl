"""
$(TYPEDSIGNATURES)

Computes the minimizer of the geometric Freidlin-Wentzell action based on the geometric
minimum action method (gMAM), using optimizers of OPtimization.jl or the original formulation
by Heymann and Vanden-Eijnden[^1]. Only Freidlin-Wentzell action has a geometric formulation.

To set an initial path different from a straight line, see the multiple dispatch method

  - `geometric_min_action_method(sys::CoupledSDEs, init::Matrix, arclength::Float64; kwargs...)`.

## Keyword arguments

  - `maxiter::Int=100`: maximum number of optimization iterations before the alogrithm stops
  - `abstol=1e-8`: absolute tolerance of action gradient to determine convergence
  - `reltol=1e-8`: relative tolerance of action gradient to determine convergence
  - `method = Adam()`: minimization algorithm (see [`Optimization.jl`](https://docs.sciml.ai/Optimization/stable/optimization_packages/optimisers/))
  - `=0.1`: step size parameter in gradient descent HeymannVandenEijnden implementation.
  - `verbose=false`: if true, print additional output
  - `show_progress=true`: if true, display a progress bar

## Optimization algorithms

The `method` keyword argument takes solver methods of the
[`Optimization.jl`](https://docs.sciml.ai/Optimization/) package; alternatively,
the option `solver = "HeymannVandenEijnden"` uses the original gMAM
algorithm[^1].

[^1]: [Heymann and Vanden-Eijnden, PRL (2008)](https://link.aps.org/doi/10.1103/PhysRevLett.100.140601)
"""
function geometric_min_action_method(sys::CoupledSDEs, x_i, x_f; N=100, kwargs...)
    path = reduce(hcat, range(x_i, x_f; length=N))
    return geometric_min_action_method(sys::CoupledSDEs, path; kwargs...)
end

"""
$(TYPEDSIGNATURES)

Runs the geometric Minimum Action Method (gMAM) to find the minimum action path (instanton) from an
initial condition `init`, given a system `sys` and total arc length `arclength`.

The initial path `init` must be a matrix of size `(D, N)`, where `D` is the dimension of the
system and `N` is the number of path points.

For more information see the main method,
[`geometric_min_action_method(sys::CoupledSDEs, x_i, x_f, arclength::Float64; kwargs...)`](@ref).
"""
function geometric_min_action_method(
    sys::CoupledSDEs,
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
    path = deepcopy(init)
    N = length(init[1, :])
    alpha = zeros(N)
    arc = range(0, 1.0; length=N)

    S(x) = geometric_action(sys, fix_ends(x, init[:, 1], init[:, end]), 1.0)

    prog = Progress(maxiter; enabled=show_progress)
    if method == "HeymannVandenEijnden"
        # error("The HeymannVandenEijnden method is broken")
        @warn("The HeymannVandenEijnden method currently does not work.")
        # for i in 1:maxiter
        #     update = heymann_vandeneijnden_step(sys, path, N; tau=tau)
        #     path .= update
        #     interpolate_path!(path, alpha, arc)
        #     next!(prog)
        # end
    else
        optf = OptimizationFunction((x, _) -> S(x), AD)
        prob = OptimizationProblem(optf, init, ())

        function callback(state, loss_val)
            interpolate_path!(state.u, alpha, arc)
            show_progress ? next!(prog) : nothing
            return false
        end

        sol = solve(
            prob,
            Optimisers.Adam();
            maxiters=maxiter,
            callback=callback,
            abstol=abstol,
            reltol=reltol,
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
function heymann_vandeneijnden_step(sys::CoupledSDEs, path, N; tau=0.1, diff_order=4)
    L = 1.0
    dx = L / (N - 1)
    update = zeros(size(path))
    lambdas, lambdas_prime = zeros(N), zeros(N)
    x_prime = path_velocity(path, 0:dx:L; order=diff_order)

    A = inv(covariance_matrix(sys))

    for i in 2:(N - 1)
        lambdas[i] = anorm(drift(sys, path[:, i]), A) / anorm(path[:, i], A)
    end
    for i in 2:(N - 1)
        lambdas_prime[i] = (lambdas[i + 1] - lambdas[i - 1]) / (2 * dx)
    end

    b(x) = drift(sys, x)
    jac(x) = jacobian(sys)(x, sys.p0)

    J = [jac(path[:, i]) for i in 2:(N - 1)]
    prod1 = [(J[i - 1] - J[i - 1]') * x_prime[:, i] for i in 2:(N - 1)]
    prod2 = [(J[i - 1]') * b(path[:, i]) for i in 2:(N - 1)]

    # Solve linear system M*x = v for each system dimension
    #! might be made faster using LinearSolve.jl special solvers
    for j in 1:size(path, 1)
        M = Matrix{Float64}(1 * I(N))
        v = zeros(N)

        # Boundary conditions
        v[1] = path[j, 1]
        v[end] = path[j, end]

        # Linear system of equations
        for i in 2:(N - 1)
            alpha = tau * lambdas[i]^2 / (dx^2)
            M[i, i] += 2 * alpha
            M[i, i - 1] = -alpha
            M[i, i + 1] = -alpha

            v[i] = (
                path[j, i] - lambdas[i] * prod1[i - 1][j] - prod2[i - 1][j] +
                lambdas[i] * lambdas_prime[i] * x_prime[j, i]
            )
        end

        # Solve and store solution
        sol = M \ v
        update[j, :] = sol
    end
    return update
end
