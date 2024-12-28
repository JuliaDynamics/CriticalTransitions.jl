"""
$(TYPEDSIGNATURES)

Runs the Minimum Action Method (MAM) to find the minimum action path (instanton) between an
initial state `x_i` and final state `x_f` in phase space.

This algorithm uses the minimizers of the
[`Optimization.jl`](https://github.com/SciML/Optimization.jl) package to minimize the
action functional (see [`fw_action`](@ref)) of a path for the given CoupledSDEs
`sys`. The path is initialized as a straight line between `x_i` and `x_f`, parameterized in
time via `N` equidistant points and total time `T`. Thus, the time step between discretized
path points is ``\\Delta t = T/N``.
To set an initial path different from a straight line, see the multiple dispatch method

  - `min_action_method(sys::CoupledSDEs, init::Matrix, T::Real; kwargs...)`.

The minimization can be performed in blocks to save intermediate results.

## Keyword arguments

  - `functional = "FW"`: type of action functional to minimize.
    Defaults to [`fw_action`](@ref), alternative: [`om_action`](@ref).
  - `maxiter = 100`: maximum number of iterations before the algorithm stops.
  - `abstol=1e-8`: absolute tolerance of action gradient to determine convergence
  - `reltol=1e-8`: relative tolerance of action gradient to determine convergence
  - `method = Adam()`: minimization algorithm (see [`Optimization.jl`](https://docs.sciml.ai/Optimization/stable/optimization_packages/optimisers/))
  - `verbose = true`: whether to print Optimization information during the run
  - `show_progress = false`: whether to print a progress bar

## Output

If `save_iterations`, returns `Optim.OptimizationResults`. Else, returns only the optimizer (path).
If `blocks > 1`, a vector of results/optimizers is returned.
"""
function min_action_method(sys::CoupledSDEs, x_i, x_f, N::Int, T::Real; kwargs...)
    init = reduce(hcat, range(x_i, x_f; length=N))
    return min_action_method(sys::CoupledSDEs, init, T; kwargs...)
end;

"""
$(TYPEDSIGNATURES)

Runs the Minimum Action Method (MAM) to find the minimum action path (instanton) from an
initial condition `init`, given a system `sys` and total path time `T`.

The initial path `init` must be a matrix of size `(D, N)`, where `D` is the dimension
of the system and `N` is the number of path points. The physical time of the path
is specified by `T`, such that the time step between consecutive path points is
``\\Delta t = T/N``.

For more information see the main method,
[`min_action_method(sys::CoupledSDEs, x_i, x_f, N::Int, T::Real; kwargs...)`](@ref).
"""
function min_action_method(
    sys::CoupledSDEs,
    init::Matrix,
    T::Real;
    functional="FW",
    maxiter::Int=100,
    abstol=nothing,
    reltol=nothing,
    method=Optimisers.Adam(),
    AD=Optimization.AutoFiniteDiff(),
    verbose=false,
    show_progress=true,
)
    arc = range(0.0, T; length=size(init, 2))
    S(x) = action(sys, fix_ends(x, init[:, 1], init[:, end]), arc, functional)

    optf = OptimizationFunction((x, _) -> S(x), AD)
    prob = OptimizationProblem(optf, init, ())

    prog = Progress(maxiter; enabled=show_progress)
    function callback(state, loss_val)
        verbose && println("Loss: $loss_val")
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
    return MaximumLikelihoodPath(sol.u, sol.objective)
end;

"""
$(TYPEDSIGNATURES)

Changes the first and last row of the matrix `x` to the vectors `x_i` and `x_f`,
respectively.
"""
function fix_ends(x::Matrix, x_i, x_f)
    m = x
    m[:, 1] = x_i
    m[:, end] = x_f
    return m
end;
