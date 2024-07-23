"""
$(TYPEDSIGNATURES)

Simulates the CoupledSDEs `sys` forward in time for a duration `T`, starting at initial condition `init`.

This function uses the [`SDEProblem`](https://diffeq.sciml.ai/stable/types/sde_types/#SciMLBase.SDEProblem) functionality of `DifferentialEquations.jl`.

## Keyword arguments
* `alg=SOSRA()`: [SDE solver](https://docs.sciml.ai/DiffEqDocs/stable/solvers/sde_solve/#sde_solve). Defaults to an  adaptive strong order 1.5 method
* `kwargs...`: keyword arguments for [`solve(SDEProblem)`](https://docs.sciml.ai/DiffEqDocs/stable/basics/common_solver_opts/#solver_options)

For more info, see [`SDEProblem`](https://diffeq.sciml.ai/stable/types/sde_types/#SciMLBase.SDEProblem).

"""
function simulate(sys::CoupledSDEs, T, init=current_state(sys);
    alg=sys.integ.alg, kwargs...
)
    prob = remake(referrenced_sciml_prob(sys); u0=init, tspan=(0, T))
    return solve(prob, alg; kwargs...)
end;

"""
$(TYPEDSIGNATURES)

Simulates the deterministic dynamics of CoupledSDEs `sys` in time for a duration `T`, starting at initial condition `init`.

This function integrates `sys.f` forward in time, using the [`ODEProblem`](https://diffeq.sciml.ai/stable/types/ode_types/#SciMLBase.ODEProblem) functionality of `DifferentialEquations.jl`. Thus, `relax` is identical to [`simulate`](@ref) when setting the noise strength `sys.σ = 0`.

## Keyword arguments
* `solver=Tsit5()`: [ODE solver](https://docs.sciml.ai/DiffEqDocs/stable/solvers/ode_solve/#ode_solve). Defaults to Tsitouras 5/4 Runge-Kutta method
* `kwargs...`: keyword arguments for [`solve(ODEProblem)`](https://docs.sciml.ai/DiffEqDocs/stable/basics/common_solver_opts/#solver_options)

For more info, see [`ODEProblem`](https://diffeq.sciml.ai/stable/types/ode_types/#SciMLBase.ODEProblem).
For stochastic integration, see [`simulate`](@ref).

"""
function relax(sys::CoupledSDEs, T, init=current_state(sys);
    alg=Tsit5(), kwargs...
)
    sde_prob = referrenced_sciml_prob(sys)
    prob = ODEProblem{isinplace(sde_prob)}(dynamic_rule(sys), init, (0, T), sys.p0)
    return solve(prob, alg; kwargs...)
end;


using CriticalTransitions

function fitzhugh_nagumo(u, p, t)
    x, y = u
    ϵ, β, α, γ, κ, I = p

    dx = (-α * x^3 + γ * x - κ * y + I) / ϵ
    dy = -β * y + x

    return SA[dx, dy]
end

using DynamicalSystems
sys = CoupledSDEs(fitzhugh_nagumo, id_func, zeros(2), [1.,3.,1.,1.,1.,0.],  σ)
ls = lyapunovspectrum(CoupledODEs(sys), 10000)