using CriticalTransitions, HDF5, Dates, StatsBase, Printf

# Settings ######################################
eps = [0.01, 0.1, 0.5, 1., 10.];
sigma = [1.2,1.2,1.2,1.2,1.2];
dt = 0.001;
tmax = 1e3;#1e5
saveat = 0.002;
#################################################

println("Running $(length(eps)) trajectory simulations for $(tmax) time units at dt=$(dt) (saving every $(saveat/dt) steps) ...")

sys(ϵ, σ) = StochSystem(fitzhugh_nagumo, [ϵ,3.,1.,1.,1.,0.], 2, σ);

fp, eigs, stab = fixedpoints(sys(1.,0.), [-2,-2], [2,2]);
R, L = fp[stab];

function run_sim(eps, sigma; dt=dt, tmax=tmax, saveat=saveat, init=L, save=true)
    tstart = now()
    sim = simulate(sys(eps,sigma), init, dt=dt, tmax=tmax, saveat=saveat)
    traj = reduce(hcat, sim.u)
    runtime = canonicalize(Dates.CompoundPeriod(Millisecond(now()-tstart)))
    println(sim.retcode, ", runtime=$(runtime)")
    
    if save
        filename = "fhn_simulation_eps=$(eps)_sigma=$(sigma)_tmax=$(tmax)"
        file = make_h5(filename, "../data/")
        attributes(file)["info"] = "$(tmax) time units simulation of fitzhugh_nagumo system for time scale parameter $(eps) at noise amplitude $(sigma)"
        write(file, "trajectory", traj)
        write(file, "initial_condition", Vector(L))
        attributes(file)["sys_info"] = sys_string(sys(eps, sigma))
        attributes(file)["run_info"] = "Run info:\n* simulation time step: dt=$(dt)\n* dataset time step: Dt=$(saveat) (saved every $(round(Int, saveat/dt)). data point)\n* total time: tmax=$(tmax)\n* solver: Euler()\n* runtime: $(runtime)"    
        close(file)
        println("Completed eps=$(eps), sigma=$(sigma), closed file $(filename)")
    else
        return traj, runtime
    end
end;

function get_density(eps, sigma; dt=dt, tmax=tmax, saveat=saveat, init=L, bins=100, urange=(-2,2), vrange=(-1.5,1.5), save=true)
    traj, runtime = run_sim(eps, sigma; dt=dt, tmax=tmax, saveat=saveat, init=L, save=false)

    hist = fit(Histogram, (traj[1,:], traj[2,:]), (-2:4/bins:2, -1.5:3/bins:1.5))
    bin_area = (4/bins)*(3/bins)
    total_counts = size(traj, 2)
    density = hist.weights' / total_counts / bin_area

    if save
        filename = @sprintf("fhn_density_eps%.2e_sigma%.2e_N%.1e", eps, sigma, total_counts)
        println(filename)
        file = make_h5(filename)#, "../data/")
        attributes(file)["info"] = "Density rho(u,v) obtained from $(tmax) time units simulation of fitzhugh_nagumo system for time scale parameter $(eps) at noise amplitude $(sigma)"
        write(file, "density", density)
        write(file, "binedges_u", Vector(hist.edges[1]))
        write(file, "binedges_v", Vector(hist.edges[2]))
        attributes(file)["sys_info"] = sys_string(sys(eps, sigma))
        attributes(file)["run_info"] = "Run info:\n* simulation time step: dt=$(dt)\n* dataset time step: Dt=$(saveat) (saved every $(round(Int, saveat/dt)). data point)\n* total time: tmax=$(tmax)\n* initial condition: $(Vector(L))\n* solver: Euler()\n* runtime: $(runtime)"
        attributes(file)["data_info"] = "Additional dataset info:\n* dimensions (1,2) = (u,v)\n* nbins=$(bins)\n* density defined as counts per total counts per bin area"
        close(file)
    else
        hist.edges[1], hist.edges[2], density
    end
end

for i in 1:length(eps)
    get_density(eps[i],sigma[i], save=true)
end;

println("Done.")