include("../StochSystem.jl")

"""
    equilib(sys::StochSystem, state::State; kwargs...)
Returns the equilibrium solution of the system `sys` for given initial condition `state`.

> Warning: This algorithm simply evolves the deterministic system forward in time until a steady-state condition is satisfied.
> Thus, the algorithm may output a false solution if it gets stuck in a quasi-equilibrium, or slowly evolving state.
> For more robust results, use `fixedpoints`.

## Keyword arguments:
* `abstol = 1e-5`: steady-state condition. Simulation ends when the rate of change (Euclidean distance in state space) of the state falls below `abstol`.
* `tmax = 1e5`: maximum simulation time before the algorithm stops even if the steady-state condition is not reached.
* `dt = 0.01`: time step of the ODE solver.
* `solver = Euler()`: ODE solver used for evolving the state.
"""
function equilib(sys::StochSystem, state::State;
    dt=0.01,
    tmax=1e5,
    abstol=1e-5,
    solver=Euler())
    
    condition(u, t, integrator) = norm(integrator.uprev-u) < abstol
    affect!(integrator) = terminate!(integrator)
    equilib_cond = DiscreteCallback(condition, affect!)

    prob = ODEProblem(sys.f, state, (0, tmax), p(sys))
    sol = solve(prob, solver; dt=dt, callback=equilib_cond, save_on=false, save_start=false)

    sol.u[1]
end;

function fixedpoints(sys::StochSystem, box)
    DynamicalSystems.fixedpoints(tocds(sys), box)
end;

function fixedpoints(sys::StochSystem, bmin::Vector, bmax::Vector)
    box = intervals_to_box(bmin, bmax)
    DynamicalSystems.fixedpoints(tocds(sys), box)
end;