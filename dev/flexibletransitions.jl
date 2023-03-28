function transition2(sys::StochSystem, x_i::State, x_f1::State, x_f2::State;
    rad_i=0.1,
    rad_f=0.1,
    dt=0.01,
    tmax=1e3,
    solver=EM(),
    progress=true,
    cut_start=true,
    rad_dims=1:length(sys.u),
    kwargs...)

    condition(u,t,integrator) = subnorm(u - x_f1; directions=rad_dims) < rad_f || subnorm(u - x_f2; directions=rad_dims) < rad_f
    affect!(integrator) = terminate!(integrator)
    cb_ball = DiscreteCallback(condition, affect!)

    sim = simulate(sys, x_i, dt=dt, tmax=tmax, solver=solver, callback=cb_ball, progress=progress, kwargs...)

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

function transitions2(sys::StochSystem, x_i::State, x_f1::State, x_f2::State, N=1;
    rad_i=0.1,
    rad_f=0.1,
    dt=0.01,
    tmax=1e3,
    Nmax=1000,
    solver=EM(),
    cut_start=true,
    rad_dims=1:length(sys.u),
    savefile=nothing,
    showprogress::Bool=true)
    """
    Generates N transition samples of sys from x_i to x_f.
    Supports multi-threading.
    rad_i:      ball radius around x_i
    rad_f:      ball radius around x_f
    cut_start:  if false, saves the whole trajectory up to the transition
    savefile:   if not nothing, saves data to a specified open .jld2 file
    """

    samples, times, idx::Vector{Int64}, r_idx::Vector{Int64} = [], [], [], []

    iterator = showprogress ? tqdm(1:Nmax) : 1:Nmax

    Threads.@threads for j ∈ iterator
        
        sim, simt, success = transition2(sys, x_i, x_f1, x_f2;
                    rad_i=rad_i, rad_f=rad_f, rad_dims=rad_dims, dt=dt, tmax=tmax,
                    solver=solver, progress=false, cut_start=cut_start)
        
        if success 
            
            if showprogress
                print("\rStatus: $(length(idx))/$(N) transitions complete.")
            end

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
        else
            push!(r_idx, j)
        end
    end

    samples, times, idx, length(r_idx)
end;