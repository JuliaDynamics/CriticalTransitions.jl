using DifferentialEquations, ProgressBars

include("StochSystem.jl")
include("noise.jl")
include("io.jl")

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
        sol
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

    sol
end;

function transition(sys::StochSystem, x_i::State, x_f::State;
    rad_i=1.0,
    rad_f=1.0,
    dt=0.01,
    tmax=1e2,
    solver=EM(),
    progress=true,
    cut_start=true)
    """
    Simulates a sample transition from x_i to x_f.
    rad_i:      ball radius around x_i
    rad_f:      ball radius around x_f
    cut_start:  if false, saves the whole trajectory up to the transition
    """

    condition(u,t,integrator) = norm(u - x_f) < rad_f
    affect!(integrator) = terminate!(integrator)
    cb_ball = DiscreteCallback(condition, affect!)

    sim = simulate(sys, x_i, dt=dt, tmax=tmax, solver=solver, callback=cb_ball, progress=progress)

    success = true
    if sim.t[end] == tmax
        success = false
    end
    
    simt = sim.t
    if cut_start
        idx = size(sim)[2]
        dist = norm(sim[:,idx] - x_i)
        while dist > rad_i
            idx -= 1
            dist = norm(sim[:,idx] - x_i)
        end
        sim = sim[:,idx:end]
        simt = simt[idx:end]
    end

    sim, simt, success
end;

function transitions(sys::StochSystem, x_i::State, x_f::State, N=1;
    rad_i=1.0,
    rad_f=1.0,
    dt=0.01,
    tmax=1e2,
    solver=EM(),
    cut_start=true,
    savefile=nothing,
    progress=true,
    verbose=true)
    """
    Generates N transition samples of sys from x_i to x_f.
    Supports multi-threading.
    rad_i:      ball radius around x_i
    rad_f:      ball radius around x_f
    cut_start:  if false, saves the whole trajectory up to the transition
    savefile:   if not nothing, saves data to a specified open .jld2 file
    """

    samples, times, idx = [], [], []

    Threads.@threads for n = tqdm(1:N)
        println("Simulation $(n): starting...")
        sim, simt, success = transition(sys, x_i, x_f;
                    rad_i=rad_i, rad_f=rad_f, dt=dt, tmax=tmax,
                    solver=solver, progress=verbose, cut_start=cut_start)
        
        if success

            println("Simulation $(n): reached competing state after $(simt[end]) time units.")
            push!(idx, n)
            
            # store or save in .jld2 file
            if savefile == nothing
                push!(samples, sim);
                push!(times, simt);
            else
                write(savefile, "paths/path "*string(n), sim)
                write(savefile, "times/times "*string(n), simt)
            end
        else
            println("WARNING: Simulation $(n) stopped before a transition occurred.")
        end
    end

    samples, times, idx
end;