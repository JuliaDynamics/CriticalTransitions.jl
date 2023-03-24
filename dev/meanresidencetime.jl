#using ProgressBars, Statistics, DifferentialEquations, DynamicalSystems, IntervalRootFinding, IntervalArithmetic, LinearAlgebra, StaticArrays, CairoMakie

# function for N of sigma

function meanresidencetime(sys::StochSystem, x_i::State, x_f::State; 
    M::Int64 = 50,
    rad_i = 0.1,
    rad_f = 0.1,
    maxsamplesize = 10000)

    ## the point of this function is to return the mean residence time given a noise strength

    meanresidencetimes = zeros(M+1); # mrt

    for mrt_position ∈ 1:M
        sample_size = mrt_position;
        times = transitions(sys, x_i, x_f, sample_size; rad_i, rad_f, showprogress = false)[2];
        finaltimes = [times[i][end] for i ∈ 1:length(times)];
        meanresidencetimes[mrt_position] = mean(finaltimes);
    end

    converged = false; mrt_position = M; sample_size = M+1;

    lastMmeans = zeros(M);

    fig = Figure(fontsize = 20)

    axs1 = Axis(fig[1,1], title = "Iterative mean transition time", xlabel = "Sample size", ylabel = "Mean transition time")
    scatter!(axs1, 1:M, meanresidencetimes[1:M], color = "blue")

    axs2 = Axis(fig[2,1], title = "Convergence condition", xlabel = "Sample size", ylabel = "Max. rel. diff.") 
    hlines!(axs2, [0.025], color = "green")

    display(fig)

    while ~converged
        
        mrt_position = (mrt_position > M) ? 1 : mrt_position+1    

        #print("$(mrt_position)")

        #meanresidencetimes = push!(meanresidencetimes, 0.)

        times = transitions(sys, x_i, x_f, sample_size; rad_i, rad_f, showprogress = false)[2];
        finaltimes = [times[i][end] for i ∈ 1:length(times)];
        meanresidencetimes[mrt_position] = mean(finaltimes);

        lastMmeans .= mean(finaltimes).-meanresidencetimes[setdiff(1:M+1, [mrt_position])];
        reldiff = abs.(lastMmeans)./abs(mean(finaltimes));

        if maximum(reldiff) < 0.025 || sample_size ≥ maxsamplesize
            converged = true
        end

        print("\rStatus: Sample size $(sample_size), maximum reldiff $(maximum(reldiff)).")

        scatter!(axs1, [sample_size], [mean(finaltimes)], color = "blue")

        scatter!(axs2, [sample_size], [maximum(reldiff)], color = "red")

        display(fig)

        sample_size += 1;

    end

    return meanresidencetimes[mrt_position]

end

# α = 8.5; ξ = 1/10; σ = 0.2;
# modtb = modtb_αξσ(α, ξ, σ); # this is a stoch system!!

# bmin = [1e-12, 1e-12]; bmax = [1,1];

# fps = fixedpoints(modtb, bmin, bmax);

# sfps = fps[1][fps[3]];

# x_i = [sfps[2,1],sfps[2,2]]; x_f = [sfps[1,1], sfps[1,2]];

# meanresidencetime(modtb, x_i, x_f)

function meantransitiontime(sys::StochSystem, x_i::State, x_f::State, sample_size::Int64;
    rad_i::Float64 = 0.1,
    rad_f::Float64 = 0.1)

    times = transitions(sys, x_i, x_f, sample_size; rad_i, rad_f, showprogress = false)[2];
    transitiontimes = [times[i][end]-times[i][1] for i ∈ 1:length(times)];

    mean(transitiontimes)

end

function meanresidencetime(sys::StochSystem, x_i::State, x_f::State, sample_size::Int64; 
    rad_i::Float64 = 0.1,
    rad_f::Float64 = 0.1,
    Nmax = 1000,
    tmax = 1e3)

    times = transitions(sys, x_i, x_f, sample_size; rad_i, rad_f, tmax, Nmax)[2];
    residencetimes = [times[i][end] for i ∈ 1:length(times)];

    mean(residencetimes)

end

# ϵ = 0.1; σ = 0.125; fhn = fhn_ϵσ(ϵ, σ); 

# bmin = [-2, -2]; bmax = [2, 2];

# fps = fixedpoints(fhn, bmin, bmax);

# sfps = fps[1][fps[3]];

# x_i = [sfps[1,1], sfps[1,2]]; x_f = [sfps[2,1], sfps[2,2]];

# println("Running for σ = $σ")

# meanresidencetime(fhn, x_i, x_f, 50; Nmax = 10000, tmax = 1e4);

ϵ = 0.1; σ = 0.1; 
fhn = fhn_ϵσ(ϵ, σ);
bmin = [-2,-2]; bmax = [2, 2];
fps = fixedpoints(fhn, bmin, bmax);
sfps = fps[1][fps[3]];
x_i = [sfps[1,1], sfps[1,2]]; x_f = [sfps[2,1], sfps[2,2]];

println("Running for σ = $σ")

tr = transitions(fhn, x_i, x_f, 10; Nmax = 5000, tmax = 1e4);