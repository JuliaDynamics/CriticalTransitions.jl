using DifferentialEquations, ProgressBars

include("StochSystem.jl")
include("noise.jl")
include("io.jl")

function simulate(sys::StochSystem, init::Vector;
    dt=0.01,
    tmax=1e3,
    solver=EM(),
    callback=nothing,
    progress=true)

    if sys.process == "WhiteGauss"
        prob = SDEProblem(sys.f, σg(sys), init, (0, tmax), p(sys), noise=gauss(sys))
        sol = solve(prob, solver; dt=dt, callback=callback, progress=progress)
        sol
    else
        ArgumentError("ERROR: Noise process not yet implemented.")
    end
end;

function relax(sys::StochSystem, init::Vector;
    dt=0.01,
    tmax=1e3,
    solver=Euler(),
    callback=nothing)
    
    prob = ODEProblem(sys.f, init, (0, tmax), p(sys))
    sol = solve(prob, solver; dt=dt, callback=callback)

    sol
end;

function transition(sys::StochSystem, x_i::Vector, x_f::Vector;
    rad_i=1.0,
    rad_f=1.0,
    dt=0.01,
    tmax=1e2,
    solver=EM(),
    progress=true,
    cut_start=true)

    condition(u,t,integrator) = norm(u - x_f) < rad_f
    affect!(integrator) = terminate!(integrator)
    cb_ball = DiscreteCallback(condition, affect!)

    sim = simulate(sys, x_i, dt=dt, tmax=tmax, solver=solver, callback=cb_ball, progress=progress)

    success = true
    if sim.t[end] == tmax
        println("WARNING: Simulation stopped before a transition occurred.")
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

function transitions(sys::StochSystem, x_i::Vector, x_f::Vector, N=1;
    rad_i=1.0,
    rad_f=1.0,
    dt=0.01,
    tmax=1e2,
    solver=EM(),
    progress=true,
    cut_start=true,
    savefile=nothing)

    samples, times = [], []

    i = 1
    Threads.@threads for n = tqdm(1:N)
        sim, simt, success = transition(sys, x_i, x_f;
                    rad_i=rad_i, rad_f=rad_f, dt=dt, tmax=tmax,
                    solver=solver, progress=false, cut_start=cut_start)
        
        if success
            
            # store or save in .jld2 file
            if savefile == nothing
                push!(samples, sim);
                push!(times, simt);
            else
                write(savefile, "paths/path "*string(i), sim)
                write(savefile, "times/times "*string(i), simt)
            end

            i += 1
        end
    end

    samples, times
end;