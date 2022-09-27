using DynamicalSystems
using DifferentialEquations
using LinearAlgebra

include("StochSystem.jl")

function equilib(sys::StochSystem, state::Vector;
    dt=0.01,
    tmax=1e3,
    abstol=1e-5,
    solver=nothing)
    
    condition(u, t, integrator) = norm(integrator.uprev-u) < abstol
    affect!(integrator) = terminate!(integrator)
    equilib_cond = DiscreteCallback(condition, affect!)

    prob = ODEProblem(sys.f, state, (0, tmax), p(sys))
    sol = solve(prob, solver; dt=dt, callback=equilib_cond, save_on=false, save_start=false)

    sol.u
end;