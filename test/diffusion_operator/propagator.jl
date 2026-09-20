using CriticalTransitions
using Test

@testset "propagate_density: default Δt=T returns endpoint" begin
    sys = CoupledSDEs((u, p, t) -> [-u[1]], [0.0]; noise_strength = 1.0)
    grid = CartesianGrid((-5.0, 5.0, 100))
    gen = DiffusionGenerator(sys, grid)
    ρ_0 = zeros(100); ρ_0[50] = 1 / grid.h[1]

    ρs, t = propagate_density(gen, 5.0, ρ_0)
    @test size(ρs) == (100, 2)
    @test t == [0.0, 5.0]
    @test ρs[:, 1] == ρ_0                                        # t=0 column
    @test sum(ρs[:, 2]) * grid.h[1] ≈ 1.0 atol = 1.0e-3            # mass preserved
end

@testset "propagate_density: trajectory with Δt" begin
    sys = CoupledSDEs((u, p, t) -> [-u[1]], [0.0]; noise_strength = 1.0)
    grid = CartesianGrid((-5.0, 5.0, 100))
    gen = DiffusionGenerator(sys, grid)
    ρ_0 = zeros(100); ρ_0[50] = 1 / grid.h[1]
    ρ_inf = stationary_distribution(gen)

    ρs, t = propagate_density(gen, 10.0, ρ_0; Δt = 2.0)
    @test t == collect(0.0:2.0:10.0)
    @test size(ρs) == (100, 6)
    @test ρs[:, 1] == ρ_0
    # Distance to invariant decreases monotonically along the trajectory.
    dists = [sum(abs.(ρs[:, i] .- ρ_inf)) * grid.h[1] for i in 1:size(ρs, 2)]
    @test issorted(dists; rev = true)
    # Final snapshot has relaxed.
    @test dists[end] < 1.0e-3
end

@testset "propagate_density: Ttr shifts the recording window" begin
    sys = CoupledSDEs((u, p, t) -> [-u[1]], [0.0]; noise_strength = 1.0)
    grid = CartesianGrid((-5.0, 5.0, 100))
    gen = DiffusionGenerator(sys, grid)
    ρ_0 = zeros(100); ρ_0[50] = 1 / grid.h[1]

    ρs, t = propagate_density(gen, 5.0, ρ_0; Δt = 1.0, Ttr = 2.0)
    @test t == collect(2.0:1.0:7.0)
    @test size(ρs) == (100, 6)
end

@testset "propagate_density: stationary density is fixed" begin
    sys = CoupledSDEs((u, p, t) -> [-u[1]], [0.0]; noise_strength = 1.0)
    grid = CartesianGrid((-5.0, 5.0, 100))
    gen = DiffusionGenerator(sys, grid)
    ρ_inf = stationary_distribution(gen)

    ρs, _ = propagate_density(gen, 5.0, ρ_inf)
    @test maximum(abs.(ρs[:, end] .- ρ_inf)) < 1.0e-3
end

@testset "propagate_density: tighter tol improves accuracy" begin
    sys = CoupledSDEs((u, p, t) -> [-u[1]], [0.0]; noise_strength = 1.0)
    grid = CartesianGrid((-5.0, 5.0, 100))
    gen = DiffusionGenerator(sys, grid)
    ρ_0 = zeros(100); ρ_0[50] = 1 / grid.h[1]

    ρs_loose, _ = propagate_density(gen, 5.0, ρ_0; tol = 1.0e-3)
    ρs_tight, _ = propagate_density(gen, 5.0, ρ_0; tol = 1.0e-12, m = 60)
    # Tight should be closer to the discrete stationary distribution than loose.
    ρ_∞ = stationary_distribution(gen)
    err_loose = sum(abs, ρs_loose[:, end] .- ρ_∞) * grid.h[1]
    err_tight = sum(abs, ρs_tight[:, end] .- ρ_∞) * grid.h[1]
    @test err_tight <= err_loose
    # Mass exactly conserved with tight tolerance.
    @test sum(ρs_tight[:, end]) * grid.h[1] ≈ 1.0 atol = 1.0e-8
end

@testset "propagate_density: input validation" begin
    sys = CoupledSDEs((u, p, t) -> [-u[1]], [0.0]; noise_strength = 1.0)
    grid = CartesianGrid((-3.0, 3.0, 30))
    gen = DiffusionGenerator(sys, grid)
    @test_throws DimensionMismatch propagate_density(gen, 1.0, zeros(20))
    @test_throws ArgumentError propagate_density(gen, -1.0, zeros(30))
    @test_throws ArgumentError propagate_density(gen, 1.0, zeros(30); Δt = 0.0)
    @test_throws ArgumentError propagate_density(gen, 1.0, zeros(30); Ttr = -1.0)
end

# 1D Ornstein-Uhlenbeck `dx = -x dt + σ dW`. Starting from a Gaussian
# `(x_0, s_0²)`, the density at time t is the Gaussian
#   μ(t)   = x_0 e^{-t}
#   var(t) = s_0² e^{-2t} + (σ²/2)(1 - e^{-2t})
# Test that propagate_density reproduces this analytical evolution.
@testset "Analytical: propagator for 1D OU is exact Gaussian" begin
    σ, x_0, s_0 = 1.0, 2.0, 0.2
    sys = CoupledSDEs((u, p, t) -> [-u[1]], [0.0]; noise_strength = σ)
    grid = CartesianGrid((-5.0, 5.0, 200))
    gen = DiffusionGenerator(sys, grid)

    ρ_0 = (1 / (s_0 * sqrt(2π))) .* exp.(-(grid.centers[1] .- x_0) .^ 2 ./ (2 * s_0^2))
    ρ_0 ./= sum(ρ_0) * grid.h[1]

    for t_test in (0.5, 1.0, 2.0, 5.0)
        ρs, _ = propagate_density(gen, t_test, ρ_0; tol = 1.0e-10, m = 50)
        ρ_t = ρs[:, end]

        # Numerical mean and variance of the propagated density.
        μ_num = sum(grid.centers[1] .* ρ_t) * grid.h[1]
        var_num = sum((grid.centers[1] .- μ_num) .^ 2 .* ρ_t) * grid.h[1]

        μ_ana = x_0 * exp(-t_test)
        var_ana = s_0^2 * exp(-2 * t_test) + (σ^2 / 2) * (1 - exp(-2 * t_test))

        @test μ_num ≈ μ_ana atol = 5.0e-3
        @test var_num ≈ var_ana atol = 5.0e-3

        # Pointwise comparison to the analytical Gaussian.
        s_ana = sqrt(var_ana)
        ρ_ana =
            (1 / (s_ana * sqrt(2π))) .*
            exp.(-(grid.centers[1] .- μ_ana) .^ 2 ./ (2 * var_ana))
        L1_err = sum(abs.(ρ_t .- ρ_ana)) * grid.h[1]
        @test L1_err < 5.0e-3
    end
end

# Propagator semigroup property: exp(t1 · F) · exp(t2 · F) = exp((t1+t2) · F)
# applied to any density. A non-trivial structural test of `propagate_density`.
@testset "Analytical: propagator semigroup property" begin
    sys = CoupledSDEs((u, p, t) -> [-u[1]], [0.0]; noise_strength = 1.0)
    grid = CartesianGrid((-6.0, 6.0, 200))
    gen = DiffusionGenerator(sys, grid)

    ρ_0 = zeros(200); ρ_0[100] = 1 / grid.h[1]

    ρs1, _ = propagate_density(gen, 1.5, ρ_0; tol = 1.0e-12)
    ρs2, _ = propagate_density(gen, 1.0, ρs1[:, end]; tol = 1.0e-12)
    ρs_direct, _ = propagate_density(gen, 2.5, ρ_0; tol = 1.0e-12)

    @test maximum(abs.(ρs2[:, end] .- ρs_direct[:, end])) < 1.0e-10
end

@testset "propagate_density: mass preservation by BC" begin
    sys = CoupledSDEs((u, p, t) -> [-u[1]], [0.0]; noise_strength = 1.0)
    grid = CartesianGrid((-3.0, 3.0, 60))

    # Reflecting and Periodic both conserve mass.
    for bc in (Reflecting(), Periodic())
        gen = DiffusionGenerator(sys, grid; bc = bc)
        ρ_0 = zeros(60); ρ_0[30] = 1 / grid.h[1]
        ρs, _ = propagate_density(gen, 5.0, ρ_0; tol = 1.0e-12)
        @test sum(ρs[:, end]) * grid.h[1] ≈ 1.0 atol = 1.0e-8
    end

    # Absorbing: mass should monotonically decay.
    gen_abs = DiffusionGenerator(sys, grid; bc = Absorbing())
    ρ_0 = zeros(60); ρ_0[30] = 1 / grid.h[1]
    ts = collect(0.0:0.5:5.0)
    ρs, _ = propagate_density(gen_abs, ts[end], ρ_0; Δt = ts[2] - ts[1], tol = 1.0e-12)
    masses = [sum(ρs[:, i]) * grid.h[1] for i in 1:size(ρs, 2)]
    @test issorted(masses; rev = true)            # monotonically decreasing
    @test masses[1] ≈ 1.0 atol = 1.0e-8           # initial
    @test masses[end] < masses[1]                # actually decayed
end
