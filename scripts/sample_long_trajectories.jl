using CriticalTransitions, HDF5, Dates

# Settings ######################################
eps = [0.01, 0.1, 0.5, 1., 10.]
sigma = [0.15,0.2,0.2,0.22,0.1]
dt = 0.005
tmax = 5e0
#################################################

println("Running $(length(eps)) trajectory simulations for $(tmax) time units at dt=$(dt) ...")

sys(ϵ, σ) = StochSystem(FitzHughNagumo, [ϵ,3.,1.,1.,1.,0.], 2, σ);

fp, eigs, stab = fixedpoints(sys(1.,0.), [-2,-2], [2,2]);
R, L = fp[stab];

sim, traj, runtimes = Dict(), Dict(), Dict(), Dict()
Threads.@threads for i in 1:length(eps)
    tstart = now()
    sim[i] = simulate(sys(eps[i],sigma[i]), L, dt=dt, tmax=tmax)
    traj[i] = reduce(hcat, sim[i].u)
    runtimes[i] = canonicalize(Dates.CompoundPeriod(Millisecond(now()-tstart)))
    println(sim[i].retcode, " $(i)/$(length(5))")
end

# Save data
file = make_h5("bruteforce_all-eps", "../data/")
attributes(file)["info"] = "5 long stochastic simulations of FitzHughNagumo system for different time scale parameter values at relatively large noise"

for i in 1:length(eps)
    create_group(file, "eps=$(eps[i])")
    loc = file["eps=$(eps[i])"]
    write(loc, "trajectory", traj[i])
    write(loc, "initial_condition", Vector(L))
    attributes(loc)["sys_info"] = sys_string(sys(eps[i], sigma[i]))
    attributes(loc)["run_info"] = "Run info:\n* time step: dt=$(dt)\n* total time: tmax=$(tmax)\n* solver: Euler()\n* runtime: $(runtimes[i])"    
end

println("Done.")