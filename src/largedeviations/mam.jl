#include("../StochSystem.jl")
#include("action.jl")

"""
    mam(sys::StochSystem, x_i::State, x_f::State, N::Int, T::Real; kwargs...)
Runs the Minimum Action Method (MAM) to find the minimum action path (instanton) between an
initial state `x_i` and final state `x_f`.

This algorithm uses the minimizers of the
[`Optim`](https://julianlsolvers.github.io/Optim.jl/stable/#) package to minimize the
Freidlin-Wentzell action functional (see [`fw_action`](@ref)) for the given StochSystem
`sys`. The path is initialized as a straight line between `x_i` and `x_f`, parameterized in
time via `N` equidistant points and total time `T`. Thus, the time step between discretized
path points is `Δt = T/N`. The minimization is performed in consecutive blocks of a fixed
number of iterations each.

## Keyword arguments
* `blocks = 10`: number of blocks
* `block_iterations = 10`: number of iterations per block
* `method = LBFGS()`: minimization algorithm (see [`Optim`](https://julianlsolvers.github.io/Optim.jl/stable/#))
* `output = "MAP"`: if "MAP", outputs the final path; if "all", outputs everything
* `showprogress = true`: whether to print a progress bar

## Output
The output can be controlled via the `output` keyword argument.

## Alternative methods
* `mam(sys::StochSystem, init::Matrix, T::Real; kwargs...)`
"""
function mam(sys::StochSystem, x_i::State, x_f::State, N::Int, T::Real;
    blocks = 10,
    block_iterations = 10,
    method = LBFGS(),
    output = "MAP",
    showprogress = true)
    
    println("=== Initializing MAM action minimizer ===")

    A = inv(sys.Σ)
    init = reduce(hcat, range(x_i, x_f, length=N))
    f(x) = fw_action(sys, fix_ends(x, x_i, x_f), range(0.0, T, length=N); cov_inv=A)
    result = Vector{Optim.OptimizationResults}(undef, blocks)
    result[1] = optimize(f, init, method, Optim.Options(iterations=block_iterations))
    
    iterator = showprogress ? tqdm(2:blocks) : 2:blocks
    for m in iterator
        result[m] = optimize(f, result[m-1].minimizer, method,
            Optim.Options(iterations=block_iterations))
    end

    if output == "MAP"
        return Optim.minimizer(result[end])
    elseif output == "all"
        return result
    end
end

"""
    mam(sys::StochSystem, init::Matrix, T::Real; kwargs...)
Alternative method for `mam`, taking an initial path `init` as input instead of the path's
start point, end point, and number of points.
"""
function mam(sys::StochSystem, init::Matrix, T::Real;
    blocks = 10,
    block_iterations = 10,
    method = LBFGS(),
    output = "MAP",
    showprogress = true)
    
    println("=== Initializing MAM action minimizer ===")

    A = inv(sys.Σ)
    f(x) = fw_action(sys, fix_ends(x, init[:,1], init[:,end]),
        range(0.0, T, length=size(init, 2)); cov_inv=A)
    result = Vector{Optim.OptimizationResults}(undef, blocks)
    result[1] = optimize(f, init, method, Optim.Options(iterations=block_iterations))
    
    iterator = showprogress ? tqdm(2:blocks) : 2:blocks
    for m in iterator
        result[m] = optimize(f, result[m-1].minimizer, method,
            Optim.Options(iterations=block_iterations))
    end

    if output == "MAP"
        return Optim.minimizer(result[end])
    elseif output == "all"
        return result
    end
end

"""
    fix_ends(x::Matrix, x_i::State, x_f::State)
Changes the first and last row of the matrix `x` to the vectors `x_i` and `x_f`,
respectively.
"""
function fix_ends(x::Matrix, x_i::State, x_f::State)
    m = x
    m[:,1] = x_i; m[:,end] = x_f
    m
end