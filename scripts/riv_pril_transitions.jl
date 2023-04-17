include("../src/CriticalTransitions.jl")
println("loading CT.jl")
using .CriticalTransitions
println("loading JLD2")
using JLD2

## SETTINGS
epsilons = [0.01, 1.0];
sigmas = [0.022,0.017,0.065,0.068];
run = [false, true, false, false]
start = [0.7667174092898188, 0.8640753583928289];
N = 50;
dt = 1e-2;
tmax = 1e3;
Nmax = 100000;
rad_i = 0.03;
rad_f = 0.1;
###########

g(u,p,t) = multiplicative_idx(u,p,t, [true, true]);
sys(ϵ, σ) = StochSystem(CriticalTransitions.rivals, [ϵ, 0.1, 0.3, 0.18, 0.1], zeros(2), σ, g, nothing, I(2), "WhiteGauss")

res = [];

for i in 1:length(epsilons)
    if run[2i-1]
        println("starting eps=$(epsilons[i]), 1-0, sigma=$(sigmas[2i-1])")
        @time _res = transitions(sys(epsilons[i], sigmas[2i-1]), start, [1.0,0.0], N;
            dt=dt, tmax=tmax, Nmax=Nmax, rad_i=rad_i, rad_f=rad_f, output_level=2, showprogress=true)
        save_object("230411_riv_prl_transitions_xy-x_eps$(epsilons[i])_sig$(round(sigmas[2i-1], digits=3)).jld", _res)
        push!(res, _res)
    end

    if run[2i]
        println("starting eps=$(epsilons[i]), 0-1, sigma=$(sigmas[2i])")
        @time _res = transitions(sys(epsilons[i], sigmas[2i]), start, [0.0,1.0], N;
            dt=dt, tmax=tmax, Nmax=Nmax, rad_i=rad_i, rad_f=rad_f, output_level=2, showprogress=true)
        save_object("230411_riv_prl_transitions_xy-y_eps$(epsilons[i])_sig$(round(sigmas[2i], digits=3)).jld", _res)
        push!(res, _res)
    end
end

println(res)
println("done")