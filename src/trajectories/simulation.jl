# """
#     simulate(sys::CoupledSDES, init::State; kwargs...)
# Simulates the StochSystem `sys` forward in time, starting at initial condition `init`.

# This function uses the [`SDEProblem`](https://diffeq.sciml.ai/stable/types/sde_types/#SciMLBase.SDEProblem) functionality of `DifferentialEquations.jl`.

# ## Keyword arguments
# * `dt=0.01`: time step of integration
# * `tmax=1e3`: total time of simulation
# * `solver=EM()`: [SDE solver](https://docs.sciml.ai/DiffEqDocs/stable/solvers/sde_solve/#sde_solve). Defaults to Euler-Mayurama
# * `callback=nothing`: callback condition
# * `progress=true`: shows a progress bar during simulation
# * `kwargs...`: keyword arguments for [`solve(ODEProblem)`](https://docs.sciml.ai/DiffEqDocs/stable/basics/common_solver_opts/#solver_options)

# For more info, see [`SDEProblem`](https://diffeq.sciml.ai/stable/types/sde_types/#SciMLBase.SDEProblem).

# """
function simulate(sys::CoupledSDEs, T, init = current_state(sys);
        alg = sys.integ.alg, kwargs...)
    prob = remake(referrenced_sciml_prob(sys), u0 = init, tspan = (0, T))
    solve(prob, alg; kwargs...)
end;

# """
#     relax(sys::StochSystem, init::State; kwargs...)
# Simulates the deterministic dynamics of StochSystem `sys` in time, starting at initial condition `init`.

# This function integrates `sys.f` forward in time, using the [`ODEProblem`](https://diffeq.sciml.ai/stable/types/ode_types/#SciMLBase.ODEProblem) functionality of `DifferentialEquations.jl`. Thus, `relax` is identical to [`simulate`](@ref) when setting the noise strength `sys.σ = 0`.

# ## Keyword arguments
# * `dt=0.01`: time step of integration
# * `tmax=1e3`: total time of simulation
# * `solver=Euler()`: [ODE solver](https://docs.sciml.ai/DiffEqDocs/stable/solvers/ode_solve/#ode_solve). Defaults to explicit forward Euler
# * `callback=nothing`: callback condition
# * `kwargs...`: keyword arguments for [`solve(ODEProblem)`](https://docs.sciml.ai/DiffEqDocs/stable/basics/common_solver_opts/#solver_options)

# For more info, see [`ODEProblem`](https://diffeq.sciml.ai/stable/types/ode_types/#SciMLBase.ODEProblem).
# For stochastic integration, see [`simulate`](@ref).

# """
function relax(sys::CoupledSDEs, T, init = current_state(sys); alg = Tsit5(), kwargs...)
    sde_prob = referrenced_sciml_prob(sys)
    prob = ODEProblem{isinplace(sde_prob)}(dynamic_rule(sys), init, (0, T), sys.p0)
    solve(prob, alg; kwargs...)
end;
