function benchmark_gmam_lqa!(SUITE)
    # Forced Van der Pol (Sec. 6.2 of Lin, Yu & Zhou 2018, quasi-potential ≈ 0.1567).
    function vdp(u, p, t)
        x, y = u
        return SA[x - x^3 / 3 + y - y^3 / 9, x + 0.9]
    end

    sys_det = CoupledODEs(vdp, [2.0, 0.0])
    pre = trajectory(sys_det, 200.0; Δt = 0.01)[1]
    u0 = collect(pre[end])
    sys_det = CoupledODEs(vdp, u0)
    T_period = 5.7966
    Nτ = 400
    tr = trajectory(sys_det, T_period - T_period / Nτ; Δt = T_period / Nτ)[1]
    pts = StateSpaceSet([SA[tr[k][1], tr[k][2]] for k in 1:Nτ])

    sys = CoupledSDEs(vdp, u0; noise_strength = 1.0)
    lc = LimitCycleFrame(pts, T_period, sys)
    G = local_quasipotential(lc)
    x_saddle = SA[-0.9, 0.6942]

    SUITE["Large deviation"]["Geometric minimal action"]["gMAM-LQA VdP (GeometricGradient)"] =
        @benchmarkable minimize_geometric_action(
        $sys, $lc, $x_saddle, GeometricGradient(; stepsize = 0.005);
        G = $G, tube_radius = 0.05, npoints = 160,
        inner_maxiters = 300, maxiters = 200, reltol = 1.0e-7,
    ) seconds = 60

    return nothing
end
