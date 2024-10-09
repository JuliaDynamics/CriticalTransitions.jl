"""
$(TYPEDSIGNATURES)

Simulates the CoupledSDEs `sys` forward in time for a duration `T`, starting at
initial condition `init`.

This function translates the [`CoupledSDEs`](@ref) type into an
[`SDEProblem`](https://diffeq.sciml.ai/stable/types/sde_types/#SciMLBase.SDEProblem)
and uses `solve` to integrate the system.

## Keyword arguments
* `alg=SOSRA()`: [SDE solver](https://docs.sciml.ai/DiffEqDocs/stable/solvers/sde_solve/#sde_solve). Defaults to an  adaptive strong order 1.5 method
* `kwargs...`: keyword arguments for [`solve(SDEProblem)`](https://docs.sciml.ai/DiffEqDocs/stable/basics/common_solver_opts/#solver_options)

For more info, see [`SDEProblem`](https://diffeq.sciml.ai/stable/types/sde_types/#SciMLBase.SDEProblem). 
See also [`trajectory`](@ref).
"""
function simulate(
    sys::CoupledSDEs, T, init=current_state(sys); alg=sys.integ.alg, kwargs...
)
    prob = remake(sys.integ.sol.prob; u0=init, tspan=(0, T))
    return solve(prob, alg; kwargs...)
end;

"""
$(TYPEDSIGNATURES)

Simulates the deterministic (noise-free) dynamics of CoupledSDEs `sys` in time for a
duration `T`, starting at initial condition `init`.

This function is equivalent to calling [`trajectory`](@ref) on the deterministic part of the
`CoupledSDEs` (with `noise_strength=0`). It works with the ODE solvers used for `CoupledODEs`.

## Keyword arguments
* `diffeq=(alg=Tsit5(), abstol = 1e-6, reltol = 1e-6)`: ODE solver settings (see [`CoupledODEs`](@ref))
* `kwargs...`: keyword arguments passed to [`trajectory`](@ref)

For more info, see [`ODEProblem`](https://diffeq.sciml.ai/stable/types/ode_types/#SciMLBase.ODEProblem).
For stochastic integration, see [`trajectory`](@ref) and [`simulate`](@ref).

"""
function relaxation(
    sys::CoupledSDEs, T, init=current_state(sys); diffeq=CoupledODEs(sys).diffeq, kwargs...
)
    return trajectory(CoupledODEs(sys), T, init; kwargs...)
end;
