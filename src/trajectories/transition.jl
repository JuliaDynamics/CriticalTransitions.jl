include("../StochSystem.jl")
include("../io/io.jl")
include("simulation.jl")

function transition(sys::StochSystem, x_i::State, x_f::State;
    rad_i=0.1,
    rad_f=0.1,
    dt=0.01,
    tmax=1e3,
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
    rad_i=0.1,
    rad_f=0.1,
    dt=0.01,
    tmax=1e3,
    Nmax=1000,
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

    samples, times, idx::Vector{Int64}, simt = [], [], [], []

    Threads.@threads for j = tqdm(1:Nmax)
        
        sim, simt, success = transition(sys, x_i, x_f;
                    rad_i=rad_i, rad_f=rad_f, dt=dt, tmax=tmax,
                    solver=solver, progress=false, cut_start=cut_start)
        
        if success
            print("\rStatus: $(length(idx))/$(N) transitions complete.")
        
            if savefile == nothing
                push!(samples, sim);
                push!(times, simt);
            else # store or save in .jld2/.h5 file
                write(savefile, "paths/path "*string(j), sim)
                write(savefile, "times/times "*string(j), simt)
            end
        
            push!(idx, j)

            if length(idx) > max(1, N - Threads.nthreads())
                break
            else
                continue
            end
        end
    end

    samples, times, idx
end;