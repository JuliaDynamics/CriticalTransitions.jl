"""
$(TYPEDSIGNATURES)

Computes the minimizer of the geometric Freidlin-Wentzell action based on the geometric
minimum action method (gMAM), using optimizers of OPtimization.jl or the original formulation
by Heymann and Vanden-Eijnden[^1]. Only Freidlin-Wentzell action has a geometric formulation.

To set an initial path different from a straight line, see the multiple dispatch method

  - `geometric_min_action_method(sys::CoupledSDEs, init::Matrix, arclength::Real; kwargs...)`.

## Keyword arguments

  - `maxiters::Int=100`: maximum number of optimization iterations before the algorithm stops
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
    maxiters::Int=100,
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

    prog = Progress(maxiters; enabled=show_progress)
    if method == "HeymannVandenEijnden"
        ws = heymann_vandeneijnden_workspace(sys, path)
        for i in 1:maxiters
            heymann_vandeneijnden_step!(ws, sys, path; tau=ϵ)
            path .= ws.update
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

        sol = solve(prob, method; maxiters, callback, abstol, reltol)
        path = sol.u
    end
    # if verbose && !converged
    #     @warn("Stopped after reaching maximum number of $(maxiters) iterations.")
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
    ws = heymann_vandeneijnden_workspace(sys, path)
    heymann_vandeneijnden_step!(ws, sys, path; tau=tau, diff_order=diff_order)
    return ws.update
end

struct HeymannVandenEijndenWorkspace{Tupdate,Tprime,Tvec,Tmat,Ttmp,Tarc,TA,Tjac}
    update::Tupdate
    x_prime::Tprime
    lambdas::Tvec
    lambdas_prime::Tvec
    prod1::Tmat
    prod2::Tmat
    drift_cache::Tmat
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
end

function heymann_vandeneijnden_workspace(sys::ContinuousTimeDynamicalSystem, path)
    N = size(path, 2)
    A = inv(covariance_matrix(sys))
    params = current_parameters(sys)
    jac_fun = jacobian(sys)
    jac = let jac_fun=jac_fun, params=params
        x -> jac_fun(x, params, 0.0)
    end

    update = similar(path)
    x_prime = zeros(size(path))
    lambdas = zeros(N)
    lambdas_prime = zeros(N)
    prod1 = zeros(size(path, 1), N - 2)
    prod2 = zeros(size(path, 1), N - 2)
    drift_cache = zeros(size(path, 1), N - 2)
    temp1 = zeros(size(path, 1))
    temp2 = zeros(size(path, 1))
    alpha_cache = zeros(N)
    dl = zeros(N - 1)
    d = ones(N)
    du = zeros(N - 1)
    v = zeros(N)
    arc = range(0, 1.0; length=N)

    return HeymannVandenEijndenWorkspace(
        update,
        x_prime,
        lambdas,
        lambdas_prime,
        prod1,
        prod2,
        drift_cache,
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
    )
end

function heymann_vandeneijnden_step!(
    ws::HeymannVandenEijndenWorkspace,
    sys::ContinuousTimeDynamicalSystem,
    path;
    tau=0.1,
    diff_order=4,
)
    N = size(path, 2)
    dx = 1.0 / (N - 1)

    path_velocity!(ws.x_prime, path, ws.arc; order=diff_order)

    @views for i in 2:(N - 1)
        velocity_norm = anorm(ws.x_prime[:, i], ws.A)
        if velocity_norm > 1e-14
            ws.drift_cache[:, i - 1] .= drift(sys, path[:, i])
            ws.lambdas[i] = anorm(ws.drift_cache[:, i - 1], ws.A) / velocity_norm
            if !isfinite(ws.lambdas[i])
                ws.lambdas[i] = 0.0
            end
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

    ws.d .= 1
    ws.dl .= 0
    ws.du .= 0
    @views for i in 2:(N - 1)
        ws.alpha_cache[i] = tau * ws.lambdas[i]^2 / (dx^2)
        if isfinite(ws.alpha_cache[i])
            ws.d[i] += 2 * ws.alpha_cache[i]
            ws.dl[i - 1] = -ws.alpha_cache[i]
            ws.du[i] = -ws.alpha_cache[i]
        else
            ws.alpha_cache[i] = 0.0
        end
    end

    T = LinearAlgebra.Tridiagonal(ws.dl, ws.d, ws.du)
    F = LinearAlgebra.lu(T)

    @inbounds for j in 1:size(path, 1)
        ws.v[1] = path[j, 1]
        ws.v[end] = path[j, end]
        for i in 2:(N - 1)
            if ws.alpha_cache[i] == 0.0
                ws.v[i] = path[j, i]
                continue
            end
            rhs =
                path[j, i] +
                tau * (
                    -ws.lambdas[i] * ws.prod1[j, i - 1] - ws.prod2[j, i - 1] +
                    ws.lambdas[i] * ws.lambdas_prime[i] * ws.x_prime[j, i]
                )
            ws.v[i] = isfinite(rhs) ? rhs : path[j, i]
        end
        ws.update[j, :] .= F \ ws.v
    end

    return ws.update
end
