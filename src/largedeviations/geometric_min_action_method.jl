"""
$(TYPEDSIGNATURES)

Computes the minimizer of the Freidlin-Wentzell action using the geometric minimum
action method (gMAM). Beta version, to be further documented.

To set an initial path different from a straight line, see the multiple dispatch method

  - `geometric_min_action_method(sys::CoupledSDEs, init::Matrix, arclength::Float64; kwargs...)`.

## Keyword arguments

  - `N = 100`: number of discretized path points
  - `maxiter = 100`: maximum number of iterations before the algorithm stops
  - `converge = 1e-5`: convergence threshold for absolute change in action
  - `method = LBFGS()`: choice of optimization algorithm (see below)
  - `tau = 0.1`: step size (used only if `method = "HeymannVandenEijnden"`)

## Optimization algorithms

The `method` keyword argument takes solver methods of the
[`Optim.jl`](https://julianlsolvers.github.io/Optim.jl/stable/#) package; alternatively,
the option `solver = "HeymannVandenEijnden"` uses the original gMAM
algorithm[^1].

[^1]: [Heymann and Vanden-Eijnden, PRL (2008)](https://link.aps.org/doi/10.1103/PhysRevLett.100.140601)
"""
function geometric_min_action_method(
        sys::CoupledSDEs, x_i, x_f, arclength = 1.0;
        N = 100,
        maxiter = 100,
        converge = 1e-5,
        method = LBFGS(),
        tau = 0.1,
        verbose::Bool = true,
        showprogress::Bool = false)
    path = reduce(hcat, range(x_i, x_f, length = N))
    geometric_min_action_method(sys::CoupledSDEs, path, arclength;
        maxiter = maxiter, converge = converge,
        method = method, tau = tau, verbose = verbose, showprogress = showprogress)
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

function geometric_min_action_method(sys::CoupledSDEs, init::Matrix, arclength = 1.0;
        maxiter::Int = 100,
        converge = 1e-5,
        method = LBFGS(),
        tau = 0.1,
        verbose::Bool = true,
        showprogress::Bool = false)
    verbose && println("=== Initializing gMAM action minimizer ===")

    A = inv(covariance_matrix(sys))
    path = init
    x_i = init[:, 1]
    x_f = init[:, end]
    N = length(init[1, :])

    S(x) = geometric_action(sys, fix_ends(x, x_i, x_f), arclength; cov_inv = A)
    paths = [path]
    action = [S(path)]

    iterator = showprogress ? tqdm(1:maxiter) : 1:maxiter
    for i in iterator
        if method == "HeymannVandenEijnden"
            error("The HeymannVandenEijnden method is broken")
            # update_path = heymann_vandeneijnden_step(sys, path, N, arclength;
            # tau = tau, cov_inv = A)
        else
            update = Optim.optimize(S, path, method, Optim.Options(iterations = 1))
            update_path = Optim.minimizer(update)
        end

        # re-interpolate
        path = interpolate_path(update_path, sys, N, arclength)
        push!(paths, path)
        push!(action, S(path))

        if abs(action[end] - action[end - 1]) < converge
            verbose && println("Converged after $(i) iterations.")
            return paths, action
            break
        end
    end
    verbose && @warn("Stopped after reaching maximum number of $(maxiter) iterations.")
    paths, action
end

function interpolate_path(path, sys, N, arclength)
    s = zeros(N)
    for j in 2:N
        s[j] = s[j - 1] + anorm(path[:, j] - path[:, j - 1], covariance_matrix(sys))
        #! anorm or norm?
    end
    s_length = s / s[end] * arclength
    interp = ParametricSpline(s_length, path, k = 3)
    return reduce(hcat, [interp(x) for x in range(0, arclength, length = N)])
end

"""
$(TYPEDSIGNATURES)

Solves eq. (6) of Ref.[^1] for an initial `path` with `N` points and arclength `L`.

## Keyword arguments

  - `tau = 0.1`: step size
  - `diff_order = 4`: order of the finite differencing along the path. Either `2` or `4`.
  - `cov_inv` = nothing: inverse of the covariance matrix `sys.Σ`. If `nothing`, it is computed.

[^1]: [Heymann and Vanden-Eijnden, PRL (2008)](https://link.aps.org/doi/10.1103/PhysRevLett.100.140601)
"""
function heymann_vandeneijnden_step(sys::CoupledSDEs, path, N, L;
        tau = 0.1,
        diff_order = 4,
        cov_inv = nothing)
    (cov_inv == nothing) ? A = inv(covariance_matrix(sys)) : A = cov_inv

    dx = L / (N - 1)
    update = zeros(size(path))
    lambdas, lambdas_prime = zeros(N), zeros(N)
    x_prime = path_velocity(path, 0:dx:L, order = diff_order)

    for i in 2:(N - 1)
        lambdas[i] = anorm(drift(sys, path[:, i]), A) / anorm(path[:, i], A)
    end
    for i in 2:(N - 1)
        lambdas_prime[i] = (lambdas[i + 1] - lambdas[i - 1]) / (2 * dx)
    end

    b(x) = drift(sys, x)

    J = [ForwardDiff.jacobian(b, path[:, i]) for i in 2:(N - 1)]
    prod1 = [(J[i - 1] - J[i - 1]') * x_prime[:, i] for i in 2:(N - 1)]
    prod2 = [(J[i - 1]') * b(path[:, i]) for i in 2:(N - 1)]

    # Solve linear system M*x = v for each system dimension
    #! might be made faster using LinearSolve.jl special solvers
    Threads.@threads for j in 1:size(path, 1)
        M = Matrix(1 * I(N))
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

            v[i] = (path[j, i] - lambdas[i] * prod1[i - 1][j] - prod2[i - 1][j] +
                    lambdas[i] * lambdas_prime[i] * x_prime[j, i])
        end

        # Solve and store solution
        sol = M \ v
        update[j, :] = sol
    end
    update
end
