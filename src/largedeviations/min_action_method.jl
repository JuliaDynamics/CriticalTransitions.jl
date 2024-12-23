"""
$(TYPEDSIGNATURES)

Runs the Minimum Action Method (MAM) to find the minimum action path (instanton) between an
initial state `x_i` and final state `x_f`.

This algorithm uses the minimizers of the
[`Optim`](https://julianlsolvers.github.io/Optim.jl/stable/#) package to minimize the
Freidlin-Wentzell action functional (see [`fw_action`](@ref)) for the given CoupledSDEs
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
  - `action_tol=1e-5`: relative tolerance of action value to determine convergence
  - `abstol=1e-8`: absolute tolerance of action gradient to determine convergence
  - `reltol=1e-8`: relative tolerance of action gradient to determine convergence
  - `blocks = 1`: number of iterative optimization blocks
  - `method = LBFGS()`: minimization algorithm (see [`Optim`](https://julianlsolvers.github.io/Optim.jl/stable/#))
  - `save_iterations = true`: whether to save Optim information
  - `verbose = true`: whether to print Optim information during the run
  - `show_progress = false`: whether to print a progress bar

## Output

If `save_iterations`, returns `Optim.OptimizationResults`. Else, returns only the optimizer (path).
If `blocks > 1`, a vector of results/optimizers is returned.
"""
function min_action_method(
    sys::CoupledSDEs,
    x_i,
    x_f,
    N::Int,
    T::Real;
    functional="FW",
    maxiter=100,
    blocks=1,
    method=LBFGS(),
    save_iterations=false,
    show_progress=true,
    verbose=false,
)
    init = reduce(hcat, range(x_i, x_f; length=N))
    return min_action_method(
        sys::CoupledSDEs,
        init,
        T;
        functional=functional,
        maxiter=maxiter,
        blocks=blocks,
        method=method,
        save_iterations=save_iterations,
        show_progress=show_progress,
        verbose=verbose,
    )
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
    maxiter=100,
    blocks=1,
    method=LBFGS(),
    action_tol=1e-5,
    abstol=1e-8,
    reltol=1e-8,
    save_iterations=false,
    verbose=false,
    show_progress=true,
)
    function f(x)
        return action(
            sys,
            fix_ends(x, init[:, 1], init[:, end]),
            range(0.0, T; length=size(init, 2)),
            functional,
        )
    end

    result = Vector{Optim.OptimizationResults}(undef, blocks)
    result[1] = Optim.optimize(
        f, init, method, Optim.Options(; iterations=Int(ceil(maxiter / blocks)))
    )
    verbose ? println(result[1]) : nothing

    prog = Progress(blocks; desc="Running optimization blocks... ", enabled=show_progress)
    if blocks > 1
        for m in 2:blocks
            result[m] = Optim.optimize(
                f,
                result[m - 1].minimizer,
                method,
                Optim.Options(;
                    iterations=maxiter, g_abstol=abstol, g_reltol=reltol, f_tol=action_tol
                ),
            )
            verbose ? println(result[m]) : nothing
        end
        if save_iterations
            (return result)
        else
            return [Optim.minimizer(result[i]) for i in 1:blocks]
        end
        next!(prog)
    else
        save_iterations ? (return result[end]) : return Optim.minimizer(result[end])
    end
end;

"""
$(TYPEDSIGNATURES)

Changes the first and last row of the matrix `x` to the vectors `x_i` and `x_f`,
respectively.
"""
function fix_ends(x::Matrix, x_i::Vector, x_f::Vector)
    m = x
    m[:, 1] = x_i
    m[:, end] = x_f
    return m
end;
