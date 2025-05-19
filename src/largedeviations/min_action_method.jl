"""
    min_action_method(sys::ContinuousTimeDynamicalSystem, x_i, x_f, T::Real; kwargs...)

Minimizes an action functional to obtain a minimum action path (instanton) between an
initial state `x_i` and final state `x_f` in phase space.

This algorithm uses the 
[`Optimization.jl`](https://github.com/SciML/Optimization.jl) package to minimize the
specified action functional (either [`fw_action`](@ref) or [`om_action`](@ref))
for the system `sys` over paths connecting `x_i` to `x_f` in time `T`.

The path is initialized as a straight line between `x_i` and `x_f`, parameterized in
time via `N` equidistant points and total time `T`. Thus, the time step between
discretized path points is ``\\Delta t = T/N``.
To set an initial path different from a straight line, see the multiple dispatch method

> `min_action_method(sys::ContinuousTimeDynamicalSystem, init::Matrix, T::Real; kwargs...)`.

Returns a [`MinimumActionPath`](@ref) object containing the optimized path and the action
value.

## Keyword arguments

  - `functional = "FW"`: type of action functional to minimize.
    Defaults to [`fw_action`](@ref), alternative: "OM" for [`om_action`](@ref).
  - `N = 100`: number of path points to use for the discretization of the path.
  - `noise_strength = nothing`: noise strength for the action functional. Specify only if `functional = "OM"`.
  - `method = Optimisers.Adam()`: minimization algorithm (see [`Optimization.jl`](https://docs.sciml.ai/Optimization/stable/optimization_packages/optimisers/))
  - `ad_type = Optimization.AutoFiniteDiff()`: type of automatic differentiation to use
    (see [`Optimization.jl`](https://docs.sciml.ai/Optimization/stable/optimization_packages/optimisers/))
  - `maxiter = 100`: maximum number of iterations before the algorithm stops.
  - `abstol=1e-8`: absolute tolerance of action gradient to determine convergence
  - `reltol=1e-8`: relative tolerance of action gradient to determine convergence
  - `verbose = true`: whether to print Optimization information during the run
  - `show_progress = false`: whether to print a progress bar
"""
function min_action_method(
    sys::ContinuousTimeDynamicalSystem, x_i, x_f, T::Real; N=100, kwargs...
)
    init = reduce(hcat, range(x_i, x_f; length=N))
    return min_action_method(sys, init, T; kwargs...)
end;

"""
    min_action_method(sys::ContinuousTimeDynamicalSystem, init::Matrix, T::Real; kwargs...)

Minimizes the specified action functional to obtain a minimum action path (instanton)
between fixed end points given a system `sys` and total path time `T`.

The initial path `init` must be a matrix of size `(D, N)`, where `D` is the dimension
of the system and `N` is the number of path points. The physical time of the path
is specified by `T`, such that the time step between consecutive path points is
``\\Delta t = T/N``.
"""
function min_action_method(
    sys::ContinuousTimeDynamicalSystem,
    init::Matrix{<:Real},
    T::Real;
    functional="FW",
    noise_strength=nothing,
    method=Optimisers.Adam(),
    ad_type=Optimization.AutoFiniteDiff(),
    maxiter::Int=100,
    abstol=nothing,
    reltol=nothing,
    verbose=false,
    show_progress=true,
)
    if sys isa CoupledSDEs
        proper_MAM_system(sys)
    end
    times = range(0.0, T; length=size(init, 2))
    S(x) = action(
        sys, fix_ends(x, init[:, 1], init[:, end]), times, functional; noise_strength
    )

    optf = OptimizationFunction((x, _) -> S(x), ad_type)
    prob = OptimizationProblem(optf, init, ())

    prog = Progress(maxiter; enabled=show_progress)
    function callback(state, loss_val)
        verbose && println("Loss: $loss_val")
        show_progress ? next!(prog) : nothing
        return false
    end

    sol = solve(
        prob, method; maxiters=maxiter, callback=callback, abstol=abstol, reltol=reltol
    )
    return MinimumActionPath(StateSpaceSet(sol.u'), sol.objective)
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

"""
    action_minimizer(sys::ContinuousTimeDynamicalSystem, x_i, x_f, T::Real; kwargs...)

Alias for [`min_action_method`](@ref).]
"""
const action_minimizer = min_action_method
