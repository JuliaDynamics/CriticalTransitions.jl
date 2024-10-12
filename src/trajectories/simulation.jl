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
For stochastic integration, see [`trajectory`](@ref).

"""
function deterministic_orbit(
    sys::CoupledSDEs, T, init=current_state(sys); diffeq=CoupledODEs(sys).diffeq, kwargs...
)
    return trajectory(CoupledODEs(sys; diffeq), T, init; kwargs...)
end;
