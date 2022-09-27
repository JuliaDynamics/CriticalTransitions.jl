using DifferentialEquations

include("StochSystem.jl")

function relax(sys::StochSystem, init::Vector;
    dt=nothing,
    tmax=1e3,
    solver=Euler(),
    callback=nothing)
    
    prob = ODEProblem(sys.f, init, (0, tmax), sys.p)
    sol = solve(prob, solver; dt=dt, callback=callback)

    sol.u, sol.t
end;

function simulate(sys::StochSystem, init::Vector;
    dt=nothing,
    tmax=1e3,
    solver=EM(),
    callback=nothing)

    prob = SDEProblem(sys.f, sys.g, init, (0, tmax), sys.p)
    sol = solve(prob, solver; dt=dt, callback=callback)

    sol.u, sol.t
end;