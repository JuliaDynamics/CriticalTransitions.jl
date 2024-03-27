using CriticalTransitions, HDF5, Dates

# Settings ######################################
eps = [0.01, 0.1, 0.5, 1., 10.]
sigma = [1.,1.,1.,1.,1.]
dt = 0.001
tmax = 2e5
saveat = 0.01
#################################################

println("Running $(length(eps)) trajectory simulations for $(tmax) time units at dt=$(dt) ...")

sys(ϵ, σ) = StochSystem(fitzhugh_nagumo, [ϵ,3.,1.,1.,1.,0.], 2, σ);

fp, eigs, stab = fixedpoints(sys(1.,0.), [-2,-2], [2,2]);
R, L = fp[stab];

file = make_h5("bruteforce_all-eps_sigma=$(sigma[1])", "../data/")
attributes(file)["info"] = "5 long stochastic simulations of fitzhugh_nagumo system for time scale parameter $(eps) at noise amplitude $(sigma)"

sim, traj, runtimes = Dict(), Dict(), Dict()
Threads.@threads for i in 1:length(eps)
    tstart = now()
    @time sim[i] = simulate(sys(eps[i],sigma[i]), L, dt=dt, tmax=tmax, saveat=saveat)
    traj[i] = reduce(hcat, sim[i].u)
    runtimes[i] = canonicalize(Dates.CompoundPeriod(Millisecond(now()-tstart)))
    println(sim[i].retcode, " $(i)/$(length(5))")
    create_group(file, "eps=$(eps[i])")
    loc = file["eps=$(eps[i])"]
    write(loc, "trajectory", traj[i])
    write(loc, "initial_condition", Vector(L))
    attributes(loc)["sys_info"] = sys_string(sys(eps[i], sigma[i]))
    attributes(loc)["run_info"] = "Run info:\n* time step: dt=$(dt)\n* total time: tmax=$(tmax)\n* solver: Euler()\n* runtime: $(runtimes[i])"    
end

close(file)
println("Done.")