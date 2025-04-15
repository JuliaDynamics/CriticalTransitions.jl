"""
    adaptive_multilevel_splitting(sys::CoupledSDEs, start, target, score_function; kwargs...)

Runs the AMS (Adaptive Multilevel Splitting) algorithm.

## Keyword arguments
- `N=2`: number of ensemble members
- `n_kill=1`: number of trajectories to discard at each iteration
- `maxiter=100`: maximum number of iterations before the algorithm stops
"""
function adaptive_multilevel_splitting(sys::CoupledSDEs, init, score_function;
    init_score = 0.1,
    target_score = 1.0,
    N=2,
    n_kill=1,
    maxiter=100,
    tmax=1e6)

    prob = remake(sys.integ.sol.prob; u0=start, tspan=(0,tmax))

    ensemble = EnsembleProblem(prob, prob_func)
    solve(ensemble)

    return TransitionEnsemble()
end

function prob_func()
end