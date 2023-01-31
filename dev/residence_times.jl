include("../src/StochSystem.jl")
include("../src/utils.jl")
include("../src/io/io.jl")
include("../src/trajectories/simulation.jl")

function residence_time(sys::StochSystem, x_i::State, x_f::State;
    rad_i=0.1,
    rad_f=0.1,
    dt=0.01,
    tmax=1e3,
    solver=EM(),
    progress=true,
    rad_dims=1:sys.dim, 
    kwargs...)

    condition(u,t,integrator) = subnorm(u - x_f; directions=rad_dims) < rad_f
    affect!(integrator) = terminate!(integrator)
    cb_ball = DiscreteCallback(condition, affect!)

    restime = simulate(sys, x_i, dt=dt, tmax=tmax, solver=solver, callback=cb_ball, progress=progress, kwargs...).t[end]

    success = true
    if restime == tmax
        success = false
    end

    restime, success

end

function residence_times(sys::StochSystem, x_i::State, x_f::State, N=1;
    rad_i=0.1,
    rad_f=0.1,
    dt=0.01,
    tmax=1e3,
    solver=EM(),
    progress=true,
    rad_dims=1:sys.dim,
    savefile = nothing,
    showprogress::Bool = true,
    kwargs...)

    times::Vector, idx::Vector{Int64}, r_idx::Vector{Int64} = [], [], []

    iterator = showprogress ? tqdm(1:Nmax) : 1:Nmax

    Threads.@threads for j ∈ iterator
        
        restime, success = residence_time(sys, x_i, x_f;
                    rad_i, rad_f, rad_dims, dt, tmax,
                    solver, progress, kwargs...)
        
        if success 
            
            if showprogress
                print("\rStatus: $(length(idx)+1)/$(N) transitions complete.")
            end

            if savefile == nothing
                push!(times, restime);
            else # store or save in .jld2 file
                write(savefile, "times/times "*string(j), restime)
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

    times, idx, length(r_idx)

end