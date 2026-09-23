using CriticalTransitions
using LinearAlgebra: eigvals
using Test

@testset "eigenmodes (1D OU)" begin
    sys = CoupledSDEs((u, p, t) -> [-u[1]], [0.0]; noise_strength = 1.0)
    grid = CartesianGrid((-5.0, 5.0, 200))
    gen = DiffusionGenerator(sys, grid)

    λ, V = eigenmodes(gen, 5)

    @test length(λ) == 5
    @test size(V) == (200, 5)
    @test issorted(real.(λ); rev = true)
    @test abs(λ[1]) < 1.0e-10                                     # trivial mode
    @test maximum(real, V[:, 1]) - minimum(real, V[:, 1]) < 1.0e-10
    @test isapprox(real(λ[2]), -1.0; atol = 1.0e-2)                 # OU spectrum
    @test isapprox(real(λ[3]), -2.0; atol = 1.0e-2)
end

# 1D pure Brownian motion on (-L, L) with reflecting (no-flux) BCs.
# The Laplacian on this domain has spectrum λ_n = -(σ²/2) (nπ/(2L))² for
# n = 0, 1, 2, ..., with eigenfunctions cos(nπ(x+L)/(2L)).
@testset "Analytical: eigenmodes of 1D BM, reflecting BC" begin
    σ, L = 1.0, 1.0
    sys = CoupledSDEs((u, p, t) -> [0.0], [0.0]; noise_strength = σ)
    grid = CartesianGrid((-L, L, 200))
    gen = DiffusionGenerator(sys, grid)

    λ, _ = eigenmodes(gen, 5)
    for n in 0:4
        expected = -(σ^2 / 2) * (n * π / (2 * L))^2
        @test isapprox(real(λ[n + 1]), expected; atol = 0.05, rtol = 0.005)
    end
    @test all(abs.(imag.(λ)) .< 1.0e-10)
end

# 1D pure Brownian motion on the periodic ring (-π, π). The Laplacian on
# the circle has eigenvalues -(σ²/2) k² for k = 0, 1, 2, ..., with the
# k = 0 mode unique (constant) and k ≥ 1 doubly degenerate (cos kx, sin kx).
@testset "Analytical: eigenmodes of 1D BM, periodic BC" begin
    σ = 1.0
    sys = CoupledSDEs((u, p, t) -> [0.0], [0.0]; noise_strength = σ)
    grid = CartesianGrid((-π, π, 200))
    gen = DiffusionGenerator(sys, grid; bc = Periodic())

    # Use DenseEigen here so the doubly-degenerate pairs match to dense
    # accuracy; the default iterative backend converges each eigenpair
    # independently and won't match them to machine precision.
    λ, _ = eigenmodes(gen, 7, DenseEigen())
    expected = [0.0, -σ^2 / 2, -σ^2 / 2, -2σ^2, -2σ^2, -4.5σ^2, -4.5σ^2]
    for n in 1:7
        @test isapprox(real(λ[n]), expected[n]; atol = 1.0e-2)
    end
    # Multiplicity 2 for k ≥ 1: pairs should match each other within
    # discretisation accuracy.
    @test isapprox(real(λ[2]), real(λ[3]); atol = 1.0e-6)
    @test isapprox(real(λ[4]), real(λ[5]); atol = 1.0e-6)
end

# The Fokker-Planck operator `Qᵀ` and the generator `Q` are matrix
# transposes, so they share the same spectrum (real / complex eigenvalues
# identical, eigenvectors swap left ↔ right).
@testset "Analytical: rate_matrix and fokker_planck_operator share spectrum" begin
    sys = CoupledSDEs((u, p, t) -> [u[1] - u[1]^3], [0.0]; noise_strength = 0.6)
    grid = CartesianGrid((-2.0, 2.0, 100))
    gen = DiffusionGenerator(sys, grid)

    Q = rate_matrix(gen)
    F = fokker_planck_operator(gen)

    λ_Q = sort(real.(eigvals(Matrix(Q))); rev = true)
    λ_F = sort(real.(eigvals(Matrix(F))); rev = true)
    @test maximum(abs.(λ_Q .- λ_F)) < 1.0e-9
end

# Eigenmodes of 2D OU `dx = -x dt + σ dW` are -(n+m) for n, m ≥ 0, since
# the operator is separable into two 1D OUs each with Hermite spectrum
# -k for k = 0, 1, 2, .... Multiplicity of eigenvalue -k in 2D is k+1
# (the number of integer pairs (n, m) with n + m = k).
# Slowest 6 (sorted descending): 0, -1 (×2), -2 (×3).
@testset "Analytical: 2D OU eigenmodes are separable sums" begin
    sys = CoupledSDEs((u, p, t) -> [-u[1], -u[2]], [0.0, 0.0]; noise_strength = 1.0)
    grid = CartesianGrid((-5.0, 5.0, 71), (-5.0, 5.0, 71))
    gen = DiffusionGenerator(sys, grid)

    # DenseEigen so degenerate eigenvalues match to dense accuracy.
    λ, _ = eigenmodes(gen, 6, DenseEigen())
    expected = [0.0, -1.0, -1.0, -2.0, -2.0, -2.0]
    for i in 1:6
        @test isapprox(real(λ[i]), expected[i]; atol = 0.05)
    end
    # Multiplicity-2 pair at -1 is degenerate.
    @test isapprox(real(λ[2]), real(λ[3]); atol = 1.0e-6)
    # Multiplicity-3 cluster at -2.
    @test isapprox(real(λ[4]), real(λ[5]); atol = 1.0e-6)
end

@testset "eigenmodes: k clamping and validation" begin
    sys = CoupledSDEs((u, p, t) -> [-u[1]], [0.0]; noise_strength = 1.0)
    grid = CartesianGrid((-3.0, 3.0, 30))
    gen = DiffusionGenerator(sys, grid)

    # k = 1 returns just the trivial mode.
    λ1, V1 = eigenmodes(gen, 1)
    @test length(λ1) == 1
    @test size(V1) == (30, 1)
    @test abs(λ1[1]) < 1.0e-10

    # k > N is clamped to N.
    λ_big, V_big = eigenmodes(gen, 1000)
    @test length(λ_big) == 30
    @test size(V_big) == (30, 30)

    # k < 1 errors.
    @test_throws ArgumentError eigenmodes(gen, 0)
    @test_throws ArgumentError eigenmodes(gen, -1)
end
