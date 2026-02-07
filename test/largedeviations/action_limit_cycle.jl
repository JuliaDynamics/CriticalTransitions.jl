using Test
using CriticalTransitions

@testset "Limit cycle: zero action along deterministic orbit" begin
    # Stuart–Landau (Hopf normal form) oscillator: stable limit cycle r = sqrt(μ)
    function stuart_landau(u, p, t)
        x, y = u
        μ, ω = p
        r2 = x^2 + y^2
        dx = (μ - r2) * x - ω * y
        dy = (μ - r2) * y + ω * x
        return SA[dx, dy]
    end

    μ = 1.0
    ω = 1.0
    σ = 0.2
    sys = CoupledSDEs(stuart_landau, zeros(2), (μ, ω); noise_strength=σ)

    # Two points on the same limit cycle (radius = 1)
    x_i = SA[1.0, 0.0]
    x_f = SA[0.0, 1.0]

    # Deterministic orbit segment from θ=0 to θ=π/2 with correct travel time
    Δt = (π / 2) / ω
    N = 2001
    time = range(0.0, Δt; length=N)

    path = zeros(2, N)
    for (i, t) in enumerate(time)
        path[:, i] = (sqrt(μ)) .* SA[cos(ω * t), sin(ω * t)]
    end

    @test isapprox(path[:, 1], x_i; atol=1e-12, rtol=0)
    @test isapprox(path[:, end], x_f; atol=1e-12, rtol=0)

    S_fw = fw_action(sys, path, time)
    @test S_fw ≤ 1e-10

    S_geo = geometric_action(sys, path)
    @test S_geo ≤ 1e-10
end

@testset "Limit cycle: gMAM minimizer has ~0 geometric action" begin
    function stuart_landau(u, p, t)
        x, y = u
        μ, ω = p
        r2 = x^2 + y^2
        dx = (μ - r2) * x - ω * y
        dy = (μ - r2) * y + ω * x
        return SA[dx, dy]
    end

    μ = 1.0
    ω = 1.0
    σ = 0.2
    sys = CoupledSDEs(stuart_landau, zeros(2), (μ, ω); noise_strength=σ)

    x_i = SA[1.0, 0.0]
    x_f = SA[0.0, 1.0]

    # Start from a straight-line initial path; minimizer should move it onto the LC segment.
    res = geometric_min_action_method(
        sys,
        x_i,
        x_f;
        points=500,
        maxiters=1000,
        optimizer=GeometricGradient(),
        show_progress=false,
        verbose=false,
    )

    S = res.action
    @test S ≤ 1e-6
end
