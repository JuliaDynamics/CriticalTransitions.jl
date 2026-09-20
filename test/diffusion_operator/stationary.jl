using CriticalTransitions
using Test

@testset "Invariant density (1D OU)" begin
    sys = CoupledSDEs((u, p, t) -> [-u[1]], [0.0]; noise_strength = 1.0)
    grid = CartesianGrid((-4.0, 4.0, 200))
    gen = DiffusionGenerator(sys, grid)
    ρ = stationary_distribution(gen)
    @test sum(ρ) * grid.h[1] ≈ 1.0 atol = 1.0e-12
    xs = collect(grid.centers[1])
    ρ_analytic = exp.(-xs .^ 2) ./ sqrt(pi)
    @test sum(abs.(ρ .- ρ_analytic)) * grid.h[1] < 5.0e-3
end

@testset "Invariant density (2D OU)" begin
    sys = CoupledSDEs((u, p, t) -> [-u[1], -u[2]], [0.0, 0.0]; noise_strength = 1.0)
    grid = CartesianGrid((-4.0, 4.0, 81), (-4.0, 4.0, 81))
    gen = DiffusionGenerator(sys, grid)
    ρ = stationary_distribution(gen)
    @test sum(ρ) * prod(grid.h) ≈ 1.0 atol = 1.0e-12
    xs = [grid.centers[1][I[1]] for I in CartesianIndices(grid.nbox)]
    ys = [grid.centers[2][I[2]] for I in CartesianIndices(grid.nbox)]
    ρ_analytic = vec(exp.(-(xs .^ 2 .+ ys .^ 2)) ./ pi)
    @test sum(abs.(ρ .- ρ_analytic)) * prod(grid.h) < 5.0e-3
end

@testset "stationary_distribution: solver kwargs forwarded to LinearSolve" begin
    using LinearSolve: KrylovJL_GMRES
    sys = CoupledSDEs((u, p, t) -> [u[1] - u[1]^3], [0.0]; noise_strength = 0.6)
    grid = CartesianGrid((-2.0, 2.0, 200))
    gen = DiffusionGenerator(sys, grid)
    @test stationary_distribution(gen, KrylovJL_GMRES(); abstol = 1.0e-14, reltol = 1.0e-14) ≈
        stationary_distribution(gen) atol = 1.0e-10
end

@testset "Analytical: variance of stationary 1D OU = σ²/2" begin
    σ = 1.0
    sys = CoupledSDEs((u, p, t) -> [-u[1]], [0.0]; noise_strength = σ)
    grid = CartesianGrid((-6.0, 6.0, 200))
    gen = DiffusionGenerator(sys, grid)
    ρ = stationary_distribution(gen)
    μ = sum(grid.centers[1] .* ρ) * grid.h[1]
    var_num = sum((grid.centers[1] .- μ) .^ 2 .* ρ) * grid.h[1]
    @test isapprox(μ, 0.0; atol = 1.0e-10)              # symmetric → mean = 0
    @test isapprox(var_num, σ^2 / 2; atol = 1.0e-3)     # var = σ²/2
end

@testset "stationary_distribution: backends agree on well-conditioned 1D OU" begin
    sys = CoupledSDEs((u, p, t) -> [-u[1]], [0.0]; noise_strength = 1.0)
    grid = CartesianGrid((-4.0, 4.0, 100))
    gen = DiffusionGenerator(sys, grid)

    ρ_default = stationary_distribution(gen)                  # nothing → LinearSolve default
    ρ_dense = stationary_distribution(gen, DenseEigen())
    ρ_kk = stationary_distribution(gen, KrylovKitSolver())

    @test ρ_default ≈ ρ_dense rtol = 1.0e-6
    @test ρ_default ≈ ρ_kk rtol = 1.0e-6
end

@testset "quasi_stationary_distribution: backends agree" begin
    sys = CoupledSDEs((u, p, t) -> [u[1] - u[1]^3], [0.0]; noise_strength = 0.4)
    grid = CartesianGrid((-2.0, 2.0, 121))
    gen = DiffusionGenerator(sys, grid)
    right_basin = x -> x[1] > 0

    ρ_kk, λ_kk = quasi_stationary_distribution(gen, right_basin, KrylovKitSolver())
    ρ_de, λ_de = quasi_stationary_distribution(gen, right_basin, DenseEigen())

    @test λ_kk ≈ λ_de rtol = 1.0e-6
    @test ρ_kk ≈ ρ_de rtol = 1.0e-5
end

@testset "quasi_stationary_distribution / eigenmodes reject LinearSolve algs" begin
    using LinearSolve: KrylovJL_GMRES
    sys = CoupledSDEs((u, p, t) -> [u[1] - u[1]^3], [0.0]; noise_strength = 0.4)
    grid = CartesianGrid((-2.0, 2.0, 60))
    gen = DiffusionGenerator(sys, grid)
    @test_throws ArgumentError quasi_stationary_distribution(gen, x -> x[1] > 0, KrylovJL_GMRES())
    @test_throws ArgumentError eigenmodes(gen, 3, KrylovJL_GMRES())
end

@testset "quasi_stationary_distribution: 1D double-well" begin
    # Symmetric double well with V'(x) = x³ - x; metastable basins are
    # the half-lines x > 0 (right basin) and x < 0 (left basin).
    sys = CoupledSDEs((u, p, t) -> [u[1] - u[1]^3], [0.0]; noise_strength = 0.3)
    grid = CartesianGrid((-2.5, 2.5, 251))
    gen = DiffusionGenerator(sys, grid)

    right_basin = x -> x[1] > 0
    left_basin = x -> x[1] < 0

    ρ_R, λ_R = quasi_stationary_distribution(gen, right_basin)
    ρ_L, λ_L = quasi_stationary_distribution(gen, left_basin)

    # Exit rates of mirror-symmetric basins must match.
    @test λ_R ≈ λ_L rtol = 1.0e-6
    @test λ_R > 0

    # QSD is zero outside its basin, positive inside.
    @test all(ρ_R[grid.centers[1] .< 0] .== 0)
    @test all(ρ_L[grid.centers[1] .> 0] .== 0)
    @test maximum(ρ_R) > 0
    @test maximum(ρ_L) > 0

    # QSD normalised to 1 over its basin.
    @test sum(ρ_R) * prod(grid.h) ≈ 1.0 atol = 1.0e-6
    @test sum(ρ_L) * prod(grid.h) ≈ 1.0 atol = 1.0e-6

    # x → -x symmetry: the reflected right QSD equals the left QSD.
    ρ_R_flipped = ρ_R[end:-1:1]
    @test ρ_R_flipped ≈ ρ_L rtol = 1.0e-5
end

@testset "quasi_stationary_distribution: empty basin throws" begin
    sys = CoupledSDEs((u, p, t) -> [-u[1]], [0.0]; noise_strength = 1.0)
    grid = CartesianGrid((-3.0, 3.0, 30))
    gen = DiffusionGenerator(sys, grid)
    @test_throws ArgumentError quasi_stationary_distribution(gen, x -> false)
end

@testset "stationary_distribution: diagnostic warning fires on metastable generator" begin
    # Maier-Stein-like with two FP basins and very small noise: the
    # discrete generator has two near-zero eigenvalues within machine
    # precision; the result depends on the pin row.
    sys = CoupledSDEs((u, p, t) -> [u[1] - u[1]^3], [0.0]; noise_strength = 0.03)
    grid = CartesianGrid((-2.5, 2.5, 401))
    gen = DiffusionGenerator(sys, grid)
    @test_logs (:warn, r"may be unreliable"i) stationary_distribution(gen; verbose = true)
    # Probe is off by default — cheap checks alone don't fire on this case.
    @test_nowarn stationary_distribution(gen)
end
