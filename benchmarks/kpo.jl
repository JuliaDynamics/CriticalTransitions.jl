using OrdinaryDiffEqLowOrderRK: Euler, RK4

function benchmark_KPO!(SUITE)
    λ = 3 / 1.21 * 2 / 295
    ω0 = 1.0
    ω = 1.0
    γ = 1 / 295
    η = 0
    α = -1

    function fu(u, v)
        return (-4 * γ * ω * u - 2 * λ * v - 4 * (ω0 - ω^2) * v - 3 * α * v * (u^2 + v^2)) /
            (8 * ω)
    end
    function fv(u, v)
        return (-4 * γ * ω * v - 2 * λ * u + 4 * (ω0 - ω^2) * u + 3 * α * u * (u^2 + v^2)) /
            (8 * ω)
    end
    stream(u, v) = Point2f(fu(u, v), fv(u, v))
    dfvdv(u, v) = (-4 * γ * ω + 6 * α * u * v) / (8 * ω)
    dfudu(u, v) = (-4 * γ * ω - 6 * α * u * v) / (8 * ω)
    dfvdu(u, v) = (-2 * λ + 4 * (ω0 - ω^2) + 9 * α * u^2 + 3 * α * v^2) / (8 * ω)
    dfudv(u, v) = (-2 * λ - 4 * (ω0 - ω^2) - 3 * α * u^2 - 9 * α * v^2) / (8 * ω)

    function H_x!(out, x, p) # in-place ℜ² → ℜ²
        @inbounds @views begin
            u = x[1, :]; v = x[2, :]
            pu = p[1, :]; pv = p[2, :]
            @. out[1, :] = pu * dfudu(u, v) + pv * dfvdu(u, v)
            @. out[2, :] = pu * dfudv(u, v) + pv * dfvdv(u, v)
        end
        return nothing
    end
    function H_p!(out, x, p) # in-place ℜ² → ℜ²
        @inbounds @views begin
            u = x[1, :]; v = x[2, :]
            pu = p[1, :]; pv = p[2, :]
            @. out[1, :] = pu + fu(u, v)
            @. out[2, :] = pv + fv(u, v)
        end
        return nothing
    end

    sys = FreidlinWentzellHamiltonian{true, 2}(H_x!, H_p!)

    Nt = 100  # number of discrete time steps
    s = collect(range(0; stop = 1, length = Nt))

    xa = [-0.0208, 0.0991]
    xb = -xa
    xsaddle = [0.0, 0.0]

    xx = @. (xb[1] - xa[1]) * s + xa[1] + 4 * s * (1 - s) * xsaddle[1]
    yy = @. (xb[2] - xa[2]) * s + xa[2] + 4 * s * (1 - s) * xsaddle[2] + 0.01 * sin(2π * s)
    x_initial = Matrix([xx yy]')

    sgmam_opt = GeometricGradient(; max_backtracks = 0, stepsize = 10.0e2)
    SUITE["Large deviation"]["Simple geometric minimal action"]["KPO"] = @benchmarkable minimize_geometric_action(
        $sys, $x_initial, $sgmam_opt; maxiters = 10_000, show_progress = false
    ) seconds = 10

    SUITE["Large deviation"]["String method"]["Kerr parametric resonator"] = @benchmarkable string_method(
        $sys,
        $x_initial;
        maxiters = 10_000,
        stepsize = 0.5,
        show_progress = false,
        integrator = Euler(),
    ) seconds = 10

    # function KPO_SA(x, p, t)
    #     u, v = x
    #     return SA[fu(u, v), fv(u, v)]
    # end
    # ds_sa = CoupledSDEs(KPO_SA, zeros(2), ())

    # SUITE["Large deviation"]["Geometric minimal action"]["KPO (Optimisers.Adam; AutoFiniteDiff)"] = @benchmarkable minimize_geometric_action(
    #     $ds_sa, $x_initial; maxiters=100, show_progress=false
    # ) seconds = 20
    return nothing
end
