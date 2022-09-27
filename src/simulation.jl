using DifferentialEquations

include("StochSystem.jl")

function simulate(sys::StochSystem, init::Vector;
    dt=nothing,
    tmax=1e3,
    solver=EM(),
    callback=nothing)

    prob = SDEProblem(sys.f, σg(sys), init, (0, tmax), p(sys))
    sol = solve(prob, solver; dt=dt, callback=callback)

    sol
end;

function relax(sys::StochSystem, init::Vector;
    dt=0.01,
    tmax=1e3,
    solver=Euler(),
    callback=nothing)
    
    prob = ODEProblem(sys.f, init, (0, tmax), p(sys))
    sol = solve(prob, solver; dt=dt, callback=callback)

    sol.u, sol.t
end;