struct TransitionPathEnsemble
    paths::Vector
    times::Vector
    success_rate::Real
    t_res::Real
    t_trans::Real
end;

function prettyprint(tpe::TransitionPathEnsemble)
    return "Transition path ensemble of $(length(tpe.times)) samples
           - sampling success rate:      $(round(tpe.success_rate, digits=3))
           - mean residence time:        $(round(tpe.t_res, digits=3))
           - mean transition time:       $(round(tpe.t_trans, digits=3))
           - normalized transition rate: $(round(tpe.t_res/tpe.t_trans, digits=1))"
end

Base.show(io::IO, tpe::TransitionPathEnsemble) = print(io, prettyprint(tpe))

"""
$(TYPEDSIGNATURES)

Generates a sample transition from point `x_i` to point `x_f`.

This function simulates `sys` in time, starting from initial condition `x_i`, until entering a `length(sys.u)`-dimensional ball of radius `rad_f` around `x_f`.

## Keyword arguments
* `rad_i=0.1`: radius of ball around `x_i`
* `rad_f=0.1`: radius of ball around `x_f`
* `tmax=1e3`: maximum time when the simulation stops even `x_f` has not been reached
* `rad_dims=1:length(sys.u)`: the directions in phase space to consider when calculating the radii
  `rad_i` and `rad_f`. Defaults to all directions. To consider only a subspace of state space,
  insert a vector of indices of the dimensions to be included.
* `cut_start=true`: if `false`, returns the whole trajectory up to the transition
* `kwargs...`: keyword arguments passed to [`CommonSolve.solve`](https://docs.sciml.ai/DiffEqDocs/stable/basics/common_solver_opts/#CommonSolve.solve-Tuple{SciMLBase.AbstractDEProblem,%20Vararg{Any}})

## Output
`[path, times, success]`
* `path` (Matrix): transition path (size [dim × N], where N is the number of time points)
* `times` (Vector): time values (since start of simulation) of the path points (size N)
* `success` (bool): if `true`, a transition occured (i.e. the ball around `x_f` has been reached), else `false`

See also [`transitions`](@ref), [`trajectory`](@ref).
"""
function transition(
    sys::CoupledSDEs,
    x_i,
    x_f;
    rad_i=0.1,
    rad_f=0.1,
    tmax=1e3,
    rad_dims=1:length(current_state(sys)),
    cut_start=true,
    kwargs...,
)
    condition(u, t, integrator) = subnorm(u - x_f; directions=rad_dims) < rad_f
    affect!(integrator) = terminate!(integrator)
    cb_ball = DiscreteCallback(condition, affect!)

    prob = remake(sys.integ.sol.prob; u0=x_i, tspan=(0, tmax))
    sim = solve(prob, sys.integ.alg; callback=cb_ball, kwargs...)
    success = sim.retcode == SciMLBase.ReturnCode.Terminated

    simt = sim.t
    if cut_start
        idx = size(sim)[2]
        dist = norm(sim[:, idx] - x_i)
        while dist > rad_i
            idx -= 1
            dist = norm(sim[:, idx] - x_i)
            idx < 1 && error(
                "Trajactory never left the initial state sphere. Increase tmax or decrease rad_i.",
            )
        end
        sim = sim[:, idx:end]
        simt = simt[idx:end]
    end

    return sim, simt, success
end;

"""
    function transitions(sys::CoupledSDEs, x_i, x_f, N=1; kwargs...)

Generates an ensemble of `N` transition samples of `sys` from point `x_i` to point `x_f`.

This function repeatedly calls the [`transition`](@ref) function to efficiently generate an ensemble of transitions. Multi-threading is enabled.

## Keyword arguments
  - `rad_i=0.1`: radius of ball around `x_i`
  - `rad_f=0.1`: radius of ball around `x_f`
  - `Nmax`: number of attempts before the algorithm stops even if less than `N` transitions occurred.
  - `tmax=1e3`: maximum time when the simulation stops even `x_f` has not been reached
  - `rad_dims=1:length(sys.u)`: the directions in phase space to consider when calculating the radii
    `rad_i` and `rad_f`. Defaults to all directions. To consider only a subspace of state space,
    insert a vector of indices of the dimensions to be included.
  - `cut_start=true`: if `false`, returns the whole trajectory up to the transition
  - `showprogress`: shows a progress bar with respect to `Nmax`
  - `kwargs...`: keyword arguments passed to [`CommonSolve.solve`](https://docs.sciml.ai/DiffEqDocs/stable/basics/common_solver_opts/#CommonSolve.solve-Tuple{SciMLBase.AbstractDEProblem,%20Vararg{Any}})

See also [`transition`](@ref).

## Output

`[samples, times, idx, N_fail]`

  - `samples` (Array of Matrices): sample paths. Each path i has size (dim × N_i), where N_i is the number of path points
  - `times` (Array of Vectors): time values (since simulation start) of path points for each path
  - `idx` (Array): list of sample indices i that produced a transition
  - `N_fail` (Int): number of samples that failed to produce a transition

> An example script using `transitions` is available [here](https://github.com/juliadynamics/CriticalTransitions.jl/blob/main/examples/sample_transitions_h5.jl).
"""
function transitions(
    sys::CoupledSDEs,
    x_i,
    x_f,
    N=1;
    rad_i=0.1,
    rad_f=0.1,
    tmax=1e3,
    Nmax=1000,
    cut_start=true,
    rad_dims=1:length(current_state(sys)),
    showprogress::Bool=false,
    kwargs...,
)

    samples, times, idx::Vector{Int64}, r_idx::Vector{Int64} = [], [], [], []

    # iterator = showprogress ? tqdm(1:Nmax) : 1:Nmax

    Threads.@threads for j in 1:Nmax
        sim, simt, success = transition(
            sys,
            x_i,
            x_f;
            rad_i=rad_i,
            rad_f=rad_f,
            rad_dims=rad_dims,
            tmax=tmax,
            cut_start=cut_start,
            kwargs...,
        )

        if success
            if showprogress
                print("\rStatus: $(length(idx)+1)/$(N) transitions complete.")
            end

            push!(samples, sim)
            push!(times, simt)
            push!(idx, j)

            if length(idx) > max(1, N - Threads.nthreads())
                break
            else
                continue
            end
        else
            push!(r_idx, j)
        end
    end

    success_rate = length(idx) / (length(r_idx) + length(idx))
    mean_res_time = sum([times[i][1] for i in 1:length(times)]) + tmax * length(r_idx)
    mean_trans_time = mean([(times[i][end] - times[i][1]) for i in 1:length(times)])

    return TransitionPathEnsemble(
        samples, times, success_rate, mean_res_time, mean_trans_time
    )
end;
