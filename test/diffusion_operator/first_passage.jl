using CriticalTransitions
using Test

@testset "mean_first_passage_time" begin
    sys_ou = CoupledSDEs((u, p, t) -> [-u[1]], [0.0]; noise_strength = 1.0)
    grid_ou = CartesianGrid((-4.0, 4.0, 200))
    gen_ou = DiffusionGenerator(sys_ou, grid_ou)
    τ = mean_first_passage_time(gen_ou, x -> abs(x[1]) > 2)
    @test all(τ .>= 0)
    target_mask = [abs(x) > 2 for x in grid_ou.centers[1]]
    @test maximum(abs.(τ[target_mask])) < 1.0e-10
    @test minimum(τ[.!target_mask]) > 0
    i0 = findfirst(x -> x >= 0, grid_ou.centers[1])
    @test τ[i0] ≈ maximum(τ) atol = 1.0e-3 * maximum(τ)
end

@testset "first_passage_variance (1D OU)" begin
    sys = CoupledSDEs((u, p, t) -> [-u[1]], [0.0]; noise_strength = 1.0)
    grid = CartesianGrid((-4.0, 4.0, 200))
    gen = DiffusionGenerator(sys, grid)
    target = x -> abs(x[1]) > 2

    var = first_passage_variance(gen, target)
    τ = mean_first_passage_time(gen, target)

    target_mask = [abs(x) > 2 for x in grid.centers[1]]
    @test all(var[.!target_mask] .>= 0)            # variance ≥ 0 on free cells
    @test maximum(abs, var[target_mask]) < 1.0e-10   # zero on the target

    # σ/μ ≈ 1 for a barrier-crossing process where the mean is dominated
    # by the rare-event timescale.
    i0 = argmin(abs.(grid.centers[1]))
    cv = sqrt(var[i0]) / τ[i0]
    @test 0.5 < cv < 1.5
end

@testset "Krylov solver via LinearSolve" begin
    using LinearSolve: KrylovJL_GMRES
    sys = CoupledSDEs((u, p, t) -> [u[1] - u[1]^3], [0.0]; noise_strength = 0.6)
    grid = CartesianGrid((-2.0, 2.0, 200))
    gen = DiffusionGenerator(sys, grid)
    B = x -> x[1] > 0.7

    alg = KrylovJL_GMRES()

    @test stationary_distribution(gen, alg) ≈ stationary_distribution(gen) atol = 1.0e-10
    @test mean_first_passage_time(gen, B; alg = alg) ≈ mean_first_passage_time(gen, B) atol = 1.0e-8
end

# Pure 1D Brownian motion `dx = σ dW` on (-L, L) with reflecting BCs.
# With target = boundary cells, the BVP `(σ²/2) T'' = -1, T(±L) = 0` has
# the parabolic exit-time profile T(x) = (L² - x²) / σ². The first-passage
# *variance* at x = 0 is Var = 2L⁴ / (3σ⁴) (closed form via the second
# moment equation `(σ²/2) T₂'' = -2T`).
@testset "Analytical: MFPT and variance for 1D BM in (-L, L)" begin
    σ, L = 1.0, 1.0
    sys = CoupledSDEs((u, p, t) -> [0.0], [0.0]; noise_strength = σ)
    grid = CartesianGrid((-L, L, 401))
    gen = DiffusionGenerator(sys, grid)

    # Effective L is the cell-center position of the boundary cell.
    L_eff = L - grid.h[1] / 2
    target = x -> abs(x[1]) > L_eff - 1.0e-12

    τ = mean_first_passage_time(gen, target)
    var = first_passage_variance(gen, target)

    # T(x) = (L_eff² - x²) / σ² over the free cells.
    for x_test in (0.0, 0.3, 0.5, 0.7)
        i = argmin(abs.(grid.centers[1] .- x_test))
        x = grid.centers[1][i]
        @test τ[i] ≈ (L_eff^2 - x^2) / σ^2 atol = 1.0e-3
    end

    # Var[τ | x = 0] = 2 L_eff⁴ / (3 σ⁴) and σ_τ/μ_τ = √(2/3).
    i0 = argmin(abs.(grid.centers[1]))
    @test var[i0] ≈ 2 * L_eff^4 / (3 * σ^4) atol = 1.0e-3
    @test sqrt(var[i0]) / τ[i0] ≈ sqrt(2 / 3) atol = 1.0e-4
end

@testset "MFPT / variance: empty target" begin
    sys = CoupledSDEs((u, p, t) -> [-u[1]], [0.0]; noise_strength = 1.0)
    grid = CartesianGrid((-2.0, 2.0, 50))
    gen = DiffusionGenerator(sys, grid)

    @test_throws ArgumentError mean_first_passage_time(gen, x -> false)
    @test_throws ArgumentError first_passage_variance(gen, x -> false)
end

@testset "Type stability" begin
    sys = CoupledSDEs((u, p, t) -> [-u[1]], [0.0]; noise_strength = 1.0)
    grid = CartesianGrid((-4.0, 4.0, 50))
    gen = DiffusionGenerator(sys, grid)
    @test (@inferred stationary_distribution(gen)) isa Vector{Float64}
    @test (@inferred mean_first_passage_time(gen, x -> abs(x[1]) > 2)) isa Vector{Float64}
end
