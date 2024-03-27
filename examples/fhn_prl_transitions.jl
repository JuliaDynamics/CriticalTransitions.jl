include("../src/CriticalTransitions.jl")
println("loading CT.jl")
using .CriticalTransitions
println("loading JLD2")
using JLD2

## SETTINGS
epsilons = [0.01, 0.1, 1.0, 10.0];
sigmas = [0.0784, 0.119, 0.1199, 0.08];
start = [sqrt(2/3), sqrt(2/27)];
final = -start;
N = 100;
dt = [2e-3, 1e-2, 1e-2, 1e-2];
tmax = 1e3;
Nmax = 100000;
rad_i = 0.03;
rad_f = 0.1;
###########

sys(eps, sig) = StochSystem(fitzhugh_nagumo, [eps, 3., 1., 1., 1., 0.], zeros(2), sig)

res = Dict()

for i in 1:length(epsilons)
    println("starting eps=$(epsilons[i])")
    @time res[epsilons[i]] = transitions(sys(epsilons[i], sigmas[i]), start, final, N;
        dt=dt[i], tmax=tmax, Nmax=Nmax, rad_i=rad_i, rad_f=rad_f, output_level=2)
    save_object("230411_fhn_prl_transitions_eps$(epsilons[i]).jld", res[epsilons[i]])
end

println("done")