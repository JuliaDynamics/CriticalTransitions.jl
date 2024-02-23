include("../noiseprocesses/stochprocess.jl")

"""
    simulate(sys::StochSystem, init::State; kwargs...)
Simulates the StochSystem `sys` forward in time, starting at initial condition `init`.

This function uses the [`SDEProblem`](https://diffeq.sciml.ai/stable/types/sde_types/#SciMLBase.SDEProblem) functionality of `DifferentialEquations.jl`.

## Keyword arguments
* `dt=0.01`: time step of integration
* `tmax=1e3`: total time of simulation
* `solver=EM()`: [SDE solver](https://docs.sciml.ai/DiffEqDocs/stable/solvers/sde_solve/#sde_solve). Defaults to Euler-Mayurama
* `callback=nothing`: callback condition
* `progress=true`: shows a progress bar during simulation
* `kwargs...`: keyword arguments for [`solve(ODEProblem)`](https://docs.sciml.ai/DiffEqDocs/stable/basics/common_solver_opts/#solver_options)

For more info, see [`SDEProblem`](https://diffeq.sciml.ai/stable/types/sde_types/#SciMLBase.SDEProblem).

> Warning: This function has only been tested for the `EM()` solver and out-of-place `SDEFunction`s.
"""
function simulate(sys::StochSystem, init::State;
    dt=0.01,
    tmax=1e3,
    solver=EM(),
    callback=nothing,
    showprogress=false,
    kwargs...)

    prob = SDEProblem(sys.f, σg(sys), init, (0, tmax), p(sys), noise=stochprocess(sys))
    solve(prob, solver; dt=dt, callback=callback, progress=showprogress, kwargs...)
end;

"""
    relax(sys::StochSystem, init::State; kwargs...)
Simulates the deterministic dynamics of StochSystem `sys` in time, starting at initial condition `init`.

This function integrates `sys.f` forward in time, using the [`ODEProblem`](https://diffeq.sciml.ai/stable/types/ode_types/#SciMLBase.ODEProblem) functionality of `DifferentialEquations.jl`. Thus, `relax` is identical to [`simulate`](@ref) when setting the noise strength `sys.σ = 0`.

## Keyword arguments
* `dt=0.01`: time step of integration
* `tmax=1e3`: total time of simulation
* `solver=Euler()`: [ODE solver](https://docs.sciml.ai/DiffEqDocs/stable/solvers/ode_solve/#ode_solve). Defaults to explicit forward Euler
* `callback=nothing`: callback condition
* `kwargs...`: keyword arguments for [`solve(ODEProblem)`](https://docs.sciml.ai/DiffEqDocs/stable/basics/common_solver_opts/#solver_options)

For more info, see [`ODEProblem`](https://diffeq.sciml.ai/stable/types/ode_types/#SciMLBase.ODEProblem).
For stochastic integration, see [`simulate`](@ref).

> Warning: This function has only been tested for the `Euler()` solver.
"""
function relax(sys::StochSystem, init::State;
    dt=0.01,
    tmax=1e3,
    solver=Euler(),
    callback=nothing,
    kwargs...)

    prob = ODEProblem(sys.f, init, (0, tmax), p(sys))
    solve(prob, solver; dt=dt, callback=callback, kwargs...)
end;
