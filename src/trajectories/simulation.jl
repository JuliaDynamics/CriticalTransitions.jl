include("../StochSystem.jl")
include("../noiseprocesses/gaussian.jl")

function simulate(sys::StochSystem, init::State;
    dt=0.01,
    tmax=1e3,
    solver=EM(),
    callback=nothing,
    progress=true)
    """
    Integrates sys forward in time
    init: initial condition
    """

    if sys.process == "WhiteGauss"
        prob = SDEProblem(sys.f, σg(sys), init, (0, tmax), p(sys), noise=gauss(sys))
        sol = solve(prob, solver; dt=dt, callback=callback, progress=progress)
    else
        ArgumentError("ERROR: Noise process not yet implemented.")
    end
end;

function relax(sys::StochSystem, init::State;
    dt=0.01,
    tmax=1e3,
    solver=Euler(),
    callback=nothing)
    """
    Integrates sys forward in time in the absence of noise (deterministic dynamics).
    init: initial condition
    """
    
    prob = ODEProblem(sys.f, init, (0, tmax), p(sys))
    sol = solve(prob, solver; dt=dt, callback=callback)
end;