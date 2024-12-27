struct TransitionPathEnsemble{SSS,T,Tstat,ES}
    paths::Vector{SSS}
    times::Vector{T}
    success_rate::Tstat
    residence_time::Tstat
    transition_time::Tstat
    sciml_ensemble::ES
end;

function prettyprint(tpe::TransitionPathEnsemble)
    return "Transition path ensemble of $(length(tpe.times)) samples
           - sampling success rate:      $(round(tpe.success_rate, digits=3))
           - mean residence time:        $(round(tpe.residence_time, digits=3))
           - mean transition time:       $(round(tpe.transition_time, digits=3))
           - normalized transition rate: $(round(tpe.residence_time/tpe.transition_time, digits=1))"
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
    radius::NTuple{2}=(0.1, 0.1),
    tmax=1e3,
    radius_directions=1:length(current_state(sys)),
    cut_start=true,
    kwargs...,
)
    rad_i, rad_f = radius
    prob, cb_ball = prepare_transition_problem(
        sys, (x_i, x_f), radius, radius_directions, tmax)

    sim = solve(prob, trajectory_algorithm(sys); callback=cb_ball, kwargs...)
    success = sim.retcode == SciMLBase.ReturnCode.Terminated

    if success && cut_start
        sim = remove_start(sim, x_i, rad_i)
    end

    return StateSpaceSet(sim.u), sim.t, success
end

function prepare_transition_problem(sys, x, radius, rad_dims, tmax)
    x_i, x_f = x
    _, rad_f = radius
    condition(u, t, integrator) = subnorm(u - x_f; directions=rad_dims) < rad_f
    affect!(integrator) = terminate!(integrator)
    cb_ball = DiscreteCallback(condition, affect!)
    prob = remake(sys.integ.sol.prob; u0=x_i, tspan=(0, tmax))
    return prob, cb_ball
end

remove_start(sol, idx) = SciMLBase.setproperties(sol; u=sol.u[idx:end], t=sol.t[idx:end])
function remove_start(sol, x_i, rad_i)
    idx = size(sol)[2]
    dist = norm(sol[:, idx] - x_i)
    while dist > rad_i
        idx -= 1
        dist = norm(sol[:, idx] - x_i)
        idx < 1 && error(
            "Trajactory never left the initial state sphere. Increase `tmax` or decrease `rad`.",
        )
    end
    return remove_start(sol, idx)
end

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
  - `show_progress`: shows a progress bar with respect to `Nmax`
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
    N::Int=1;
    radius::NTuple{2}=(0.1, 0.1),
    tmax=1e3,
    Nmax=100,
    cut_start=true,
    rad_dims=1:length(current_state(sys)),
    show_progress::Bool=true,
    EnsembleAlg=EnsembleThreads()::SciMLBase.BasicEnsembleAlgorithm,
    kwargs...,
)
    prob, cb_ball = prepare_transition_problem(sys, (x_i, x_f), radius, rad_dims, tmax)

    tries = 0
    success = 0
    function output_func(sol, i)
        rerun = sol.retcode != SciMLBase.ReturnCode.Terminated && i < Nmax
        tries += 1
        !rerun && (success += 1)
        if !rerun && cut_start
            sol = remove_start(sol, x[1], radius[1])
        end
        return (sol, rerun)
    end
    ensemble = EnsembleProblem(prob; output_func=output_func)
    sim = solve(
        ensemble,
        trajectory_algorithm(sys),
        EnsembleAlg;
        callback=cb_ball,
        trajectories=N,
        kwargs...,
    )

    success_rate = success / tries
    mean_res_time = mean([sol.t[1] for sol in sim])
    mean_trans_time = mean([(sol.t[end] - sol.t[1]) for sol in sim])

    samples = [StateSpaceSet(sol.u) for sol in sim]
    times = [sol.t for sol in sim]

    return TransitionPathEnsemble(
        samples, times, success_rate, mean_res_time, mean_trans_time, sim
    )
end;
