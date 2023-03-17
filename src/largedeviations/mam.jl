"""
    mam(sys::StochSystem, x_i::State, x_f::State, N::Int, T::Real; kwargs...)
Runs the Minimum Action Method (MAM) to find the minimum action path (instanton) between an
initial state `x_i` and final state `x_f`.

This algorithm uses the minimizers of the
[`Optim`](https://julianlsolvers.github.io/Optim.jl/stable/#) package to minimize the
Freidlin-Wentzell action functional (see [`fw_action`](@ref)) for the given StochSystem
`sys`. The path is initialized as a straight line between `x_i` and `x_f`, parameterized in
time via `N` equidistant points and total time `T`. Thus, the time step between discretized
path points is ``\\Δt = T/N``.
To set an initial path different from a straight line, see the multiple dispatch method

* `mam(sys::StochSystem, init::Matrix, T::Real; kwargs...)`.

The minimization can be performed in blocks to save intermediate results.

## Keyword arguments
* `functional = "FW"`: type of action functional to minimize.  
   Defaults to [`fw_action`](@ref), alternative: [`om_action`](@ref).
* `maxiter = 100`: maximum number of iterations before the algorithm stops.
* `blocks = 1`: number of iterative optimization blocks
* `method = LBFGS()`: minimization algorithm (see [`Optim`](https://julianlsolvers.github.io/Optim.jl/stable/#))
* `save_info = true`: whether to save Optim information
* `verbose = true`: whether to print Optim information during the run
* `showprogress = false`: whether to print a progress bar
* `kwargs...`: any keyword arguments from `Optim.Options` (see [docs](http://julianlsolvers.github.io/Optim.jl/stable/#user/config/))

## Output
If `save_info`, returns `Optim.OptimizationResults`. Else, returns only the optimizer (path).
If `blocks > 1`, a vector of results/optimizers is returned.
"""
function mam(sys::StochSystem, x_i::State, x_f::State, N::Int, T::Real;
    functional="FW",
    maxiter = 100,
    blocks = 1,
    method = LBFGS(),
    save_info = true,
    showprogress = false,
    verbose = true,
    kwargs...)
    
    init = reduce(hcat, range(x_i, x_f, length=N))
    mam(sys::StochSystem, init, T;
        functional=functional, maxiter=maxiter, blocks=blocks, method=method,
        save_info=save_info, showprogress=showprogress, verbose=verbose, kwargs...)
end;

"""
    mam(sys::StochSystem, init::Matrix, T::Real; kwargs...)
Runs the Minimum Action Method (MAM) to find the minimum action path (instanton) from an
initial condition `init`, given a system `sys` and total path time `T`.

The initial path `init` must be a matrix of size `(D, N)`, where `D` is the dimension
`sys.dim` of the system and `N` is the number of path points. The physical time of the path
is specified by `T`, such that the time step between consecutive path points is
``\\Deltat = T/N``.

For more information see the main method,
[`mam(sys::StochSystem, x_i::State, x_f::State, N::Int, T::Real; kwargs...)`](@ref).
"""
function mam(sys::StochSystem, init::Matrix, T::Real;
    functional="FW",
    maxiter = 100,
    blocks = 1,
    method = LBFGS(),
    save_info = true,
    showprogress = false,
    verbose = true,
    kwargs...)
    
    println("=== Initializing MAM action minimizer ===")

    A = inv(sys.Σ)
    f(x) = action(sys, fix_ends(x, init[:,1], init[:,end]),
        range(0.0, T, length=size(init, 2)), functional; cov_inv=A)
    
    result = Vector{Optim.OptimizationResults}(undef, blocks)
    result[1] = Optim.optimize(f, init, method,
        Optim.Options(iterations=Int(ceil(maxiter/blocks)), kwargs...))
    verbose ? println(result[1]) : nothing
    
    if blocks > 1
        iterator = showprogress ? tqdm(2:blocks) : 2:blocks
        for m in iterator
            result[m] = Optim.optimize(f, result[m-1].minimizer, method,
                Optim.Options(Int(ceil(maxiter/blocks)), kwargs...))
            verbose ? println(result[m]) : nothing
        end
        save_info ? (return result) : return [Optim.minimizer(result[i]) for i=1:blocks]
    else
        save_info ? (return result[end]) : return Optim.minimizer(result[end])
    end
end;

"""
    fix_ends(x::Matrix, x_i::State, x_f::State)
Changes the first and last row of the matrix `x` to the vectors `x_i` and `x_f`,
respectively.
"""
function fix_ends(x::Matrix, x_i::State, x_f::State)
    m = x
    m[:,1] = x_i; m[:,end] = x_f
    m
end;