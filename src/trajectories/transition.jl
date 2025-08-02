"""
$(TYPEDSIGNATURES)

Generates a sample transition from point `x_i` to point `x_f`.

This function simulates `sys` in time, starting from initial condition `x_i`,
until entering a ball of given radius around `x_f`. If a seed was given to `sys` the solver
is initialized with this seed. To change the seed you can pass a new seed to the `seed` keyword.

## Keyword arguments
* `radii=(0.1, 0.1)`: radius of the ball around `x_i` and `x_f`, respectively
* `tmax=1e3`: maximum time until the simulation stops even `x_f` has not been reached
* `radius_directions=1:length(sys.u)`: the directions in phase space to consider when calculating the radii
  `rad_i` and `rad_f`. Defaults to all directions. To consider only a subspace of state space,
  insert a vector of indices of the dimensions to be included.
* `cut_start=true`: if `false`, returns the whole trajectory up to the transition
* `kwargs...`: keyword arguments passed to [`CommonSolve.solve`](https://docs.sciml.ai/DiffEqDocs/stable/basics/common_solver_opts)

## Output
`[path, times, success]`
* `path` (Matrix): transition path (size [dim × N], where N is the number of time points)
* `times` (Vector): time values (since start of simulation) of the path points (size N)
* `success` (bool): if `true`, a transition occurred (i.e. the ball around `x_f` has been reached), else `false`

See also [`transitions`](@ref), [`trajectory`](@ref).
"""
function transition(
    sys::CoupledSDEs,
    x_i,
    x_f;
    radii::NTuple{2}=(0.1, 0.1),
    tmax=1e3,
    radius_directions=1:length(current_state(sys)),
    cut_start=true,
    kwargs...,
)
    rad_i, rad_f = radii
    prob, cb_ball = prepare_transition_problem(
        sys, (x_i, x_f), radii, radius_directions, tmax
    )

    sim = solve(prob, solver(sys); callback=cb_ball, kwargs...)
    success = sim.retcode == SciMLBase.ReturnCode.Terminated

    if success && cut_start
        sim = remove_start(sim, x_i, rad_i)
    end

    return StateSpaceSet(sim.u), sim.t, success
end

function prepare_transition_problem(sys, x, radii, radius_directions, tmax)
    x_i, x_f = x
    _, rad_f = radii
    condition(u, t, integrator) = subnorm(u - x_f; directions=radius_directions) < rad_f
    affect!(integrator) = terminate!(integrator)
    cb_ball = DiscreteCallback(condition, affect!)
    prob = remake(sys.integ.sol.prob; u0=x_i, tspan=(0, tmax))
    return prob, cb_ball
end

function remove_start(sol, idx)
    ConstructionBase.setproperties(sol; u=sol.u[idx:end], t=sol.t[idx:end])
end
function remove_start(sol, x_i, rad_i)
    idx = size(sol)[2]
    dist = norm(sol[:, idx] - x_i)
    while dist > rad_i
        idx -= 1
        dist = norm(sol[:, idx] - x_i)
        idx < 1 && error(
            "Trajectory never left the initial state sphere. Increase `tmax` or decrease `rad`.",
        )
    end
    return remove_start(sol, idx)
end

"""
$(TYPEDSIGNATURES)

Generates an ensemble of `N` transition samples of `sys` from point `x_i` to point `x_f`.
The transitions is by default simulated using threading. To sample the transitions in serial,
GPU or Distributed enverionment, pass the desired
[`SciMLBase.EnsembleAlgorithm`](https://docs.sciml.ai/DiffEqDocs/stable/features/ensemble/)
to the EnsembleAlg algorithm.

## Keyword arguments
  - `radii=(0.1, 0.1)`: radius of the ball around `x_i` and `x_f`, respectively
  - `Nmax`: number of attempts before the algorithm stops even if less than `N` transitions occurred.
  - `tmax=1e3`: maximum time when the simulation stops even `x_f` has not been reached
  - `radius_directions=1:length(sys.u)`: the directions in phase space to consider when calculating the radii
    `rad_i` and `rad_f`. Defaults to all directions. To consider only a subspace of state space,
    insert a vector of indices of the dimensions to be included.
  - `cut_start=true`: if `false`, returns the whole trajectory up to the transition
  - `show_progress=true`: shows a progress bar with respect to `Nmax`
  - `kwargs...`: keyword arguments passed to
    [`CommonSolve.solve`](https://docs.sciml.ai/DiffEqDocs/stable/basics/common_solver_opts)

See also [`transition`](@ref).

Returns a [`TransitionEnsemble`](@ref) object.

"""
function transitions(
    sys::CoupledSDEs,
    x_i,
    x_f,
    N::Int=1;
    radii::NTuple{2}=(0.1, 0.1),
    tmax=1e3,
    Nmax=100,
    cut_start=true,
    radius_directions=1:length(current_state(sys)),
    show_progress::Bool=true,
    EnsembleAlg=EnsembleThreads()::SciMLBase.EnsembleAlgorithm,
    kwargs...,
)
    prob, cb_ball = prepare_transition_problem(
        sys, (x_i, x_f), radii, radius_directions, tmax
    )

    tries = 0
    success = 0
    function output_func(sol, i)
        rerun = sol.retcode != SciMLBase.ReturnCode.Terminated && i < Nmax
        tries += 1
        !rerun && (success += 1)
        if !rerun && cut_start
            sol = remove_start(sol, x_i, radii[1])
        end
        return (sol, rerun)
    end

    seed = sys.integ.sol.prob.seed
    function prob_func(prob, i, repeat)
        return remake(prob; seed=rand(Random.MersenneTwister(seed + i + repeat), UInt32))
    end

    ensemble = SciMLBase.EnsembleProblem(prob; output_func=output_func, prob_func=prob_func)
    sim = solve(
        ensemble, solver(sys), EnsembleAlg; callback=cb_ball, trajectories=N, kwargs...
    )

    return TransitionEnsemble(sim, success / tries)
end;
