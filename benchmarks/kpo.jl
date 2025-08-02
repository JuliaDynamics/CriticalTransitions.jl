function benchmark_KPO!(SUITE)
    λ = 3 / 1.21 * 2 / 295
    ω0 = 1.000
    ω = 1.000
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

    function H_x(x, p) # ℜ² → ℜ²
        u, v = eachrow(x)
        pu, pv = eachrow(p)

        H_u = @. pu * dfudu(u, v) + pv * dfvdu(u, v)
        H_v = @. pu * dfudv(u, v) + pv * dfvdv(u, v)
        return Matrix([H_u H_v]')
    end
    function H_p(x, p) # ℜ² → ℜ²
        u, v = eachrow(x)
        pu, pv = eachrow(p)

        H_pu = @. pu + fu(u, v)
        H_pv = @. pv + fv(u, v)
        return Matrix([H_pu H_pv]')
    end

    sys = SgmamSystem{false,2}(H_x, H_p)

    Nt = 100  # number of discrete time steps
    s = collect(range(0; stop=1, length=Nt))

    xa = [-0.0208, 0.0991]
    xb = -xa
    xsaddle = [0.0, 0.0]

    xx = @. (xb[1] - xa[1]) * s + xa[1] + 4 * s * (1 - s) * xsaddle[1]
    yy = @. (xb[2] - xa[2]) * s + xa[2] + 4 * s * (1 - s) * xsaddle[2] + 0.01 * sin(2π * s)
    x_initial = Matrix([xx yy]')

    SUITE["Large deviation"]["Simple geometric minimal action"]["KPO"] = @benchmarkable sgmam(
        $sys, $x_initial; iterations=10_000, ϵ=10e2, show_progress=false
    ) seconds = 10

    SUITE["Large deviation"]["String method"]["Kerr parametric resonator"] = @benchmarkable string_method(
        $sys, $x_initial; iterations=10_000, ϵ=10e2, show_progress=false
    ) seconds = 10

    # function KPO_SA(x, p, t)
    #     u, v = x
    #     return SA[fu(u, v), fv(u, v)]
    # end
    # ds_sa = CoupledSDEs(KPO_SA, zeros(2), ())

    # SUITE["Large deviation"]["Geometric minimal action"]["KPO (Optimisers.Adam; AutoFiniteDiff)"] = @benchmarkable geometric_min_action_method(
    #     $ds_sa, $x_initial; maxiter=100, show_progress=false
    # ) seconds = 20
    return nothing
end
