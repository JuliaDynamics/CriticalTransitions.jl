function benchmark_maierstein!(SUITE)
    function meier_stein(u, p, t) # out-of-place
        x, y = u
        dx = x - x^3 - 10 * x * y^2
        dy = -(1 + x^2) * y
        return SA[dx, dy]
    end
    σ = 0.25
    sys = CoupledSDEs(meier_stein, zeros(2); noise_strength=σ)

    xx = range(-1.0, 1.0; length=30)
    yy = 0.3 .* (-xx .^ 2 .+ 1)
    init = Matrix([xx yy]')

    x_i = init[:, 1]
    x_f = init[:, end]

    SUITE["Large deviation"]["Geometric minimal action"]["Maier-Stein (Optimisers.Adam; AutoFiniteDiff)"] = @benchmarkable geometric_min_action_method(
        $sys, $init; maxiters=1000, show_progress=false, optimizer=Optimisers.Adam()
    ) seconds = 10
    SUITE["Large deviation"]["Geometric minimal action"]["Maier-Stein (HeymannVandenEijnden)"] = @benchmarkable geometric_min_action_method(
        $sys, $init; maxiters=1000, show_progress=false, optimizer=GeometricGradient()
    ) seconds = 10
    return nothing
end
