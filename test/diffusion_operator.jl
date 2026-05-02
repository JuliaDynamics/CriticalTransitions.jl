using Test
using LinearAlgebra: diag
using SparseArrays: rowvals, nonzeros, nzrange, SparseMatrixCSC, sparse

using CriticalTransitions

# =====================================================================
# Helpers
# =====================================================================

function _offdiag_extrema(A)
    rv = rowvals(A)
    nz = nonzeros(A)
    omin, omax = Inf, -Inf
    for col in 1:size(A, 2), p in nzrange(A, col)
        row = rv[p]
        row != col || continue
        omin = min(omin, nz[p])
        omax = max(omax, nz[p])
    end
    return omin, omax
end

# =====================================================================
# DiffusionGenerator: shape and sign conventions
# =====================================================================

@testset "DiffusionGenerator (1D OU): rate-matrix structure" begin
    sys = CoupledSDEs((u, p, t) -> [-u[1]], [0.0]; noise_strength=1.0)
    grid = CartesianGrid((-4.0, 4.0, 200))
    gen = DiffusionGenerator(sys, grid)
    @test gen isa DiffusionGenerator
    Q = gen.Q
    @test Q isa SparseMatrixCSC{Float64,Int}
    @test size(Q) == (200, 200)
    @test maximum(abs, vec(sum(Q; dims=2))) < 1e-12
    omin, omax = _offdiag_extrema(Q)
    @test omin >= 0           # rate matrix: off-diagonals = transition rates ≥ 0
    @test isfinite(omax)
    @test all(diag(Q) .<= 0)  # rate matrix: diagonal = -escape rate ≤ 0
end

@testset "rate_matrix and m_matrix accessors" begin
    sys = CoupledSDEs((u, p, t) -> [-u[1]], [0.0]; noise_strength=1.0)
    grid = CartesianGrid((-4.0, 4.0, 200))
    gen = DiffusionGenerator(sys, grid)
    Q = rate_matrix(gen)
    L = m_matrix(gen)
    @test Q === gen.Q                              # alias for the field
    @test maximum(abs, L + Q) < 1e-12              # M-matrix is the negation
    omin, omax = _offdiag_extrema(L)
    @test omax <= 0           # M-matrix: off-diagonals ≤ 0
    @test all(diag(L) .>= 0)  # M-matrix: diagonal ≥ 0
end

@testset "fokker_planck_operator" begin
    sys = CoupledSDEs((u, p, t) -> [-u[1]], [0.0]; noise_strength=1.0)
    grid = CartesianGrid((-4.0, 4.0, 50))
    gen = DiffusionGenerator(sys, grid)
    F = fokker_planck_operator(gen)
    @test F == sparse(transpose(gen.Q))
    # Stationary density solves F ρ = 0 (FP nullspace).
    ρ = stationary_distribution(gen)
    @test maximum(abs, F * ρ) < 1e-10
end

# =====================================================================
# Stationary density
# =====================================================================

@testset "Invariant density (1D OU)" begin
    sys = CoupledSDEs((u, p, t) -> [-u[1]], [0.0]; noise_strength=1.0)
    grid = CartesianGrid((-4.0, 4.0, 200))
    gen = DiffusionGenerator(sys, grid)
    ρ = stationary_distribution(gen)
    @test sum(ρ) * grid.h[1] ≈ 1.0 atol = 1e-12
    xs = collect(grid.centers[1])
    ρ_analytic = exp.(-xs .^ 2) ./ sqrt(pi)
    @test sum(abs.(ρ .- ρ_analytic)) * grid.h[1] < 5e-3
end

@testset "Invariant density (2D OU)" begin
    sys = CoupledSDEs((u, p, t) -> [-u[1], -u[2]], [0.0, 0.0]; noise_strength=1.0)
    grid = CartesianGrid((-4.0, 4.0, 81), (-4.0, 4.0, 81))
    gen = DiffusionGenerator(sys, grid)
    ρ = stationary_distribution(gen)
    @test sum(ρ) * prod(grid.h) ≈ 1.0 atol = 1e-12
    xs = [grid.centers[1][I[1]] for I in CartesianIndices(grid.nbox)]
    ys = [grid.centers[2][I[2]] for I in CartesianIndices(grid.nbox)]
    ρ_analytic = vec(exp.(-(xs .^ 2 .+ ys .^ 2)) ./ pi)
    @test sum(abs.(ρ .- ρ_analytic)) * prod(grid.h) < 5e-3
end

# =====================================================================
# Committors and MFPT (Layer 3)
# =====================================================================

@testset "Committor BCs and reversibility (1D double-well)" begin
    sys = CoupledSDEs((u, p, t) -> [u[1] - u[1]^3], [0.0]; noise_strength=0.7)
    grid = CartesianGrid((-2.0, 2.0, 200))
    gen = DiffusionGenerator(sys, grid)
    A = x -> x[1] < -0.7
    B = x -> x[1] > 0.7
    A_mask = [A((x,)) for x in grid.centers[1]]
    B_mask = [B((x,)) for x in grid.centers[1]]
    qp = forward_committor(gen, A, B)
    qm = backward_committor(gen, A, B)
    @test all(qp[A_mask] .== 0)
    @test all(qp[B_mask] .== 1)
    @test all(qm[A_mask] .== 1)
    @test all(qm[B_mask] .== 0)
    @test extrema(qp) == (0.0, 1.0)
    @test extrema(qm) == (0.0, 1.0)
    # Reversible (gradient) system: q⁻ = 1 - q⁺
    @test maximum(abs.(qm .- (1 .- qp))) < 1e-6
end

@testset "Layer-3 analyses" begin
    sys = CoupledSDEs((u, p, t) -> [u[1] - u[1]^3], [0.0]; noise_strength=0.6)
    grid = CartesianGrid((-2.0, 2.0, 200))
    gen = DiffusionGenerator(sys, grid)

    @testset "forward_committor in [0, 1]" begin
        qp = forward_committor(gen, x -> x[1] < -0.7, x -> x[1] > 0.7)
        @test extrema(qp) == (0.0, 1.0)
    end

    @testset "backward_committor with explicit reverse generator" begin
        A = x -> x[1] < -0.7
        B = x -> x[1] > 0.7
        qm_explicit = backward_committor(gen, A, B; reverse=gen)
        qm_adjoint = backward_committor(gen, A, B)
        @test maximum(abs.(qm_explicit .- qm_adjoint)) < 1e-6
    end

    @testset "mean_first_passage_time" begin
        sys_ou = CoupledSDEs((u, p, t) -> [-u[1]], [0.0]; noise_strength=1.0)
        grid_ou = CartesianGrid((-4.0, 4.0, 200))
        gen_ou = DiffusionGenerator(sys_ou, grid_ou)
        τ = mean_first_passage_time(gen_ou, x -> abs(x[1]) > 2)
        @test all(τ .>= 0)
        target_mask = [abs(x) > 2 for x in grid_ou.centers[1]]
        @test maximum(abs.(τ[target_mask])) < 1e-10
        @test minimum(τ[.!target_mask]) > 0
        i0 = findfirst(x -> x >= 0, grid_ou.centers[1])
        @test τ[i0] ≈ maximum(τ) atol = 1e-3 * maximum(τ)
    end
end

# =====================================================================
# Spectral analysis
# =====================================================================

@testset "first_passage_variance (1D OU)" begin
    sys = CoupledSDEs((u, p, t) -> [-u[1]], [0.0]; noise_strength=1.0)
    grid = CartesianGrid((-4.0, 4.0, 200))
    gen = DiffusionGenerator(sys, grid)
    target = x -> abs(x[1]) > 2

    var = first_passage_variance(gen, target)
    τ = mean_first_passage_time(gen, target)

    target_mask = [abs(x) > 2 for x in grid.centers[1]]
    @test all(var[.!target_mask] .>= 0)            # variance ≥ 0 on free cells
    @test maximum(abs, var[target_mask]) < 1e-10   # zero on the target

    # σ/μ ≈ 1 for a barrier-crossing process where the mean is dominated
    # by the rare-event timescale.
    i0 = argmin(abs.(grid.centers[1]))
    cv = sqrt(var[i0]) / τ[i0]
    @test 0.5 < cv < 1.5
end

@testset "propagate_density: default Δt=T returns endpoint" begin
    sys = CoupledSDEs((u, p, t) -> [-u[1]], [0.0]; noise_strength=1.0)
    grid = CartesianGrid((-5.0, 5.0, 100))
    gen = DiffusionGenerator(sys, grid)
    ρ_0 = zeros(100); ρ_0[50] = 1 / grid.h[1]

    ρs, t = propagate_density(gen, 5.0, ρ_0)
    @test size(ρs) == (100, 2)
    @test t == [0.0, 5.0]
    @test ρs[:, 1] == ρ_0                                        # t=0 column
    @test sum(ρs[:, 2]) * grid.h[1] ≈ 1.0 atol = 1e-3            # mass preserved
end

@testset "propagate_density: trajectory with Δt" begin
    sys = CoupledSDEs((u, p, t) -> [-u[1]], [0.0]; noise_strength=1.0)
    grid = CartesianGrid((-5.0, 5.0, 100))
    gen = DiffusionGenerator(sys, grid)
    ρ_0 = zeros(100); ρ_0[50] = 1 / grid.h[1]
    ρ_inf = stationary_distribution(gen)

    ρs, t = propagate_density(gen, 10.0, ρ_0; Δt=2.0)
    @test t == collect(0.0:2.0:10.0)
    @test size(ρs) == (100, 6)
    @test ρs[:, 1] == ρ_0
    # Distance to invariant decreases monotonically along the trajectory.
    dists = [sum(abs.(ρs[:, i] .- ρ_inf)) * grid.h[1] for i in 1:size(ρs, 2)]
    @test issorted(dists; rev=true)
    # Final snapshot has relaxed.
    @test dists[end] < 1e-3
end

@testset "propagate_density: Ttr shifts the recording window" begin
    sys = CoupledSDEs((u, p, t) -> [-u[1]], [0.0]; noise_strength=1.0)
    grid = CartesianGrid((-5.0, 5.0, 100))
    gen = DiffusionGenerator(sys, grid)
    ρ_0 = zeros(100); ρ_0[50] = 1 / grid.h[1]

    ρs, t = propagate_density(gen, 5.0, ρ_0; Δt=1.0, Ttr=2.0)
    @test t == collect(2.0:1.0:7.0)
    @test size(ρs) == (100, 6)
end

@testset "propagate_density: stationary density is fixed" begin
    sys = CoupledSDEs((u, p, t) -> [-u[1]], [0.0]; noise_strength=1.0)
    grid = CartesianGrid((-5.0, 5.0, 100))
    gen = DiffusionGenerator(sys, grid)
    ρ_inf = stationary_distribution(gen)

    ρs, _ = propagate_density(gen, 5.0, ρ_inf)
    @test maximum(abs.(ρs[:, end] .- ρ_inf)) < 1e-3
end

@testset "propagate_density: tighter tol improves accuracy" begin
    sys = CoupledSDEs((u, p, t) -> [-u[1]], [0.0]; noise_strength=1.0)
    grid = CartesianGrid((-5.0, 5.0, 100))
    gen = DiffusionGenerator(sys, grid)
    ρ_0 = zeros(100); ρ_0[50] = 1 / grid.h[1]

    ρs_loose, _ = propagate_density(gen, 5.0, ρ_0; tol=1e-3)
    ρs_tight, _ = propagate_density(gen, 5.0, ρ_0; tol=1e-12, m=60)
    @test sum(abs.(ρs_loose[:, end] .- ρs_tight[:, end])) * grid.h[1] >= 0
    # Mass exactly conserved with tight tolerance.
    @test sum(ρs_tight[:, end]) * grid.h[1] ≈ 1.0 atol = 1e-8
end

@testset "propagate_density: input validation" begin
    sys = CoupledSDEs((u, p, t) -> [-u[1]], [0.0]; noise_strength=1.0)
    grid = CartesianGrid((-3.0, 3.0, 30))
    gen = DiffusionGenerator(sys, grid)
    @test_throws DimensionMismatch propagate_density(gen, 1.0, zeros(20))
    @test_throws ArgumentError propagate_density(gen, -1.0, zeros(30))
    @test_throws ArgumentError propagate_density(gen, 1.0, zeros(30); Δt=0.0)
    @test_throws ArgumentError propagate_density(gen, 1.0, zeros(30); Ttr=-1.0)
end

@testset "eigenmodes (1D OU)" begin
    sys = CoupledSDEs((u, p, t) -> [-u[1]], [0.0]; noise_strength=1.0)
    grid = CartesianGrid((-5.0, 5.0, 200))
    gen = DiffusionGenerator(sys, grid)

    λ, V = eigenmodes(gen, 5)

    @test length(λ) == 5
    @test size(V) == (200, 5)
    @test issorted(real.(λ); rev=true)
    @test abs(λ[1]) < 1e-10                                     # trivial mode
    @test maximum(real, V[:, 1]) - minimum(real, V[:, 1]) < 1e-10
    @test isapprox(real(λ[2]), -1.0; atol=1e-2)                 # OU spectrum
    @test isapprox(real(λ[3]), -2.0; atol=1e-2)
end

# =====================================================================
# Boundary conditions
# =====================================================================

@testset "BC: Reflecting default" begin
    sys = CoupledSDEs((u, p, t) -> [-u[1]], [0.0]; noise_strength=1.0)
    grid = CartesianGrid((-4.0, 4.0, 50))
    gen = DiffusionGenerator(sys, grid)
    @test gen.bc == (Reflecting(),)
    @test gen isa DiffusionGenerator{1,Tuple{Reflecting}}
    @test maximum(abs, vec(sum(gen.Q; dims=2))) < 1e-12
end

@testset "BC: Periodic — 1D ring with pure diffusion" begin
    sys = CoupledSDEs((u, p, t) -> [0.0], [0.0]; noise_strength=1.0)
    grid = CartesianGrid((-pi, pi, 60))
    gen = DiffusionGenerator(sys, grid; bc=Periodic())
    @test gen.bc == (Periodic(),)
    @test gen isa DiffusionGenerator{1,Tuple{Periodic}}
    @test maximum(abs, vec(sum(gen.Q; dims=2))) < 1e-12
    ρ = stationary_distribution(gen)
    @test maximum(ρ) - minimum(ρ) < 1e-12
end

@testset "BC: Absorbing leaks mass on boundary" begin
    sys = CoupledSDEs((u, p, t) -> [-u[1]], [0.0]; noise_strength=1.0)
    grid = CartesianGrid((-3.0, 3.0, 60))
    gen = DiffusionGenerator(sys, grid; bc=Absorbing())
    @test gen.bc == (Absorbing(),)
    rsums = vec(sum(gen.Q; dims=2))
    @test maximum(abs, rsums[2:59]) < 1e-10
    @test rsums[1] < -1e-3
    @test rsums[end] < -1e-3
end

@testset "BC: per-axis tuple" begin
    sys = CoupledSDEs((u, p, t) -> [-u[1], 0.0], [0.0, 0.0]; noise_strength=1.0)
    grid = CartesianGrid((-2.0, 2.0, 21), (-pi, pi, 24))
    gen = DiffusionGenerator(sys, grid; bc=(Reflecting(), Periodic()))
    @test gen.bc == (Reflecting(), Periodic())
    @test gen isa DiffusionGenerator{2,Tuple{Reflecting,Periodic}}
    @test maximum(abs, vec(sum(gen.Q; dims=2))) < 1e-12
end

@testset "BC: validation errors" begin
    sys = CoupledSDEs((u, p, t) -> [-u[1]], [0.0]; noise_strength=1.0)
    grid = CartesianGrid((-3.0, 3.0, 20))
    # Symbols no longer accepted (catches typos at the type system).
    @test_throws ArgumentError DiffusionGenerator(sys, grid; bc=:reflecting)
    # Wrong tuple length.
    @test_throws ArgumentError DiffusionGenerator(
        sys, grid; bc=(Reflecting(), Periodic())
    )
end

@testset "BC: stationary_distribution rejects Absorbing" begin
    sys = CoupledSDEs((u, p, t) -> [-u[1]], [0.0]; noise_strength=1.0)
    grid = CartesianGrid((-3.0, 3.0, 30))
    gen_abs = DiffusionGenerator(sys, grid; bc=Absorbing())
    @test_throws ArgumentError stationary_distribution(gen_abs)
    # backward_committor (default adjoint route) calls stationary_distribution
    # internally, so it inherits the rejection.
    @test_throws ArgumentError backward_committor(gen_abs, x -> x[1] < -1, x -> x[1] > 1)
end

# =====================================================================
# LinearSolve Krylov backend
# =====================================================================

@testset "Krylov solver via LinearSolve" begin
    using LinearSolve: KrylovJL_GMRES
    sys = CoupledSDEs((u, p, t) -> [u[1] - u[1]^3], [0.0]; noise_strength=0.6)
    grid = CartesianGrid((-2.0, 2.0, 200))
    gen = DiffusionGenerator(sys, grid)
    A = x -> x[1] < -0.7
    B = x -> x[1] > 0.7

    alg = KrylovJL_GMRES()

    @test stationary_distribution(gen; alg=alg) ≈ stationary_distribution(gen) atol = 1e-10
    @test forward_committor(gen, A, B; alg=alg) ≈ forward_committor(gen, A, B) atol = 1e-10
    @test backward_committor(gen, A, B; alg=alg) ≈ backward_committor(gen, A, B) atol = 1e-10
    @test mean_first_passage_time(gen, B; alg=alg) ≈ mean_first_passage_time(gen, B) atol = 1e-8
end

# =====================================================================
# Analytical benchmarks: closed-form physical examples
# =====================================================================

# Pure 1D Brownian motion `dx = σ dW` on (-L, L) with reflecting BCs.
# With target = boundary cells, the BVP `(σ²/2) T'' = -1, T(±L) = 0` has
# the parabolic exit-time profile T(x) = (L² - x²) / σ². The first-passage
# *variance* at x = 0 is Var = 2L⁴ / (3σ⁴) (closed form via the second
# moment equation `(σ²/2) T₂'' = -2T`).
@testset "Analytical: MFPT and variance for 1D BM in (-L, L)" begin
    σ, L = 1.0, 1.0
    sys = CoupledSDEs((u, p, t) -> [0.0], [0.0]; noise_strength=σ)
    grid = CartesianGrid((-L, L, 401))
    gen = DiffusionGenerator(sys, grid)

    # Effective L is the cell-center position of the boundary cell.
    L_eff = L - grid.h[1] / 2
    target = x -> abs(x[1]) > L_eff - 1e-12

    τ = mean_first_passage_time(gen, target)
    var = first_passage_variance(gen, target)

    # T(x) = (L_eff² - x²) / σ² over the free cells.
    for x_test in (0.0, 0.3, 0.5, 0.7)
        i = argmin(abs.(grid.centers[1] .- x_test))
        x = grid.centers[1][i]
        @test τ[i] ≈ (L_eff^2 - x^2) / σ^2 atol = 1e-3
    end

    # Var[τ | x = 0] = 2 L_eff⁴ / (3 σ⁴) and σ_τ/μ_τ = √(2/3).
    i0 = argmin(abs.(grid.centers[1]))
    @test var[i0] ≈ 2 * L_eff^4 / (3 * σ^4) atol = 1e-3
    @test sqrt(var[i0]) / τ[i0] ≈ sqrt(2 / 3) atol = 1e-4
end

# 1D pure Brownian motion on (-L, L) with reflecting (no-flux) BCs.
# The Laplacian on this domain has spectrum λ_n = -(σ²/2) (nπ/(2L))² for
# n = 0, 1, 2, ..., with eigenfunctions cos(nπ(x+L)/(2L)).
@testset "Analytical: eigenmodes of 1D BM, reflecting BC" begin
    σ, L = 1.0, 1.0
    sys = CoupledSDEs((u, p, t) -> [0.0], [0.0]; noise_strength=σ)
    grid = CartesianGrid((-L, L, 200))
    gen = DiffusionGenerator(sys, grid)

    λ, _ = eigenmodes(gen, 5)
    for n in 0:4
        expected = -(σ^2 / 2) * (n * π / (2 * L))^2
        @test isapprox(real(λ[n + 1]), expected; atol=0.05, rtol=0.005)
    end
    @test all(abs.(imag.(λ)) .< 1e-10)
end

# 1D pure Brownian motion on the periodic ring (-π, π). The Laplacian on
# the circle has eigenvalues -(σ²/2) k² for k = 0, 1, 2, ..., with the
# k = 0 mode unique (constant) and k ≥ 1 doubly degenerate (cos kx, sin kx).
@testset "Analytical: eigenmodes of 1D BM, periodic BC" begin
    σ = 1.0
    sys = CoupledSDEs((u, p, t) -> [0.0], [0.0]; noise_strength=σ)
    grid = CartesianGrid((-π, π, 200))
    gen = DiffusionGenerator(sys, grid; bc=Periodic())

    λ, _ = eigenmodes(gen, 7)
    expected = [0.0, -σ^2 / 2, -σ^2 / 2, -2σ^2, -2σ^2, -4.5σ^2, -4.5σ^2]
    for n in 1:7
        @test isapprox(real(λ[n]), expected[n]; atol=1e-2)
    end
    # Multiplicity 2 for k ≥ 1: pairs should match each other within
    # discretisation accuracy.
    @test isapprox(real(λ[2]), real(λ[3]); atol=1e-6)
    @test isapprox(real(λ[4]), real(λ[5]); atol=1e-6)
end

# Forward committor for pure 1D BM `dx = σ dW` between A = {x_min} and
# B = {x_max} (point sets). The BVP `(σ²/2) q'' = 0` with q(x_A) = 0,
# q(x_B) = 1 gives the linear committor q(x) = (x - x_A) / (x_B - x_A).
@testset "Analytical: linear committor for 1D BM" begin
    σ, L = 1.0, 1.0
    sys = CoupledSDEs((u, p, t) -> [0.0], [0.0]; noise_strength=σ)
    grid = CartesianGrid((-L, L, 401))
    gen = DiffusionGenerator(sys, grid)

    x_A = grid.centers[1][1]
    x_B = grid.centers[1][end]
    A = x -> x[1] < x_A + 1e-12
    B = x -> x[1] > x_B - 1e-12

    qp = forward_committor(gen, A, B)
    for i in (50, 100, 200, 300, 350)
        x = grid.centers[1][i]
        @test qp[i] ≈ (x - x_A) / (x_B - x_A) atol = 1e-10
    end
end

# 1D Ornstein-Uhlenbeck `dx = -x dt + σ dW`. Starting from a Gaussian
# `(x_0, s_0²)`, the density at time t is the Gaussian
#   μ(t)   = x_0 e^{-t}
#   var(t) = s_0² e^{-2t} + (σ²/2)(1 - e^{-2t})
# Test that propagate_density reproduces this analytical evolution.
@testset "Analytical: propagator for 1D OU is exact Gaussian" begin
    σ, x_0, s_0 = 1.0, 2.0, 0.2
    sys = CoupledSDEs((u, p, t) -> [-u[1]], [0.0]; noise_strength=σ)
    grid = CartesianGrid((-5.0, 5.0, 200))
    gen = DiffusionGenerator(sys, grid)

    ρ_0 = (1 / (s_0 * sqrt(2π))) .* exp.(-(grid.centers[1] .- x_0) .^ 2 ./ (2 * s_0^2))
    ρ_0 ./= sum(ρ_0) * grid.h[1]

    for t_test in (0.5, 1.0, 2.0, 5.0)
        ρs, _ = propagate_density(gen, t_test, ρ_0; tol=1e-10, m=50)
        ρ_t = ρs[:, end]

        # Numerical mean and variance of the propagated density.
        μ_num = sum(grid.centers[1] .* ρ_t) * grid.h[1]
        var_num = sum((grid.centers[1] .- μ_num) .^ 2 .* ρ_t) * grid.h[1]

        μ_ana = x_0 * exp(-t_test)
        var_ana = s_0^2 * exp(-2 * t_test) + (σ^2 / 2) * (1 - exp(-2 * t_test))

        @test μ_num ≈ μ_ana atol = 5e-3
        @test var_num ≈ var_ana atol = 5e-3

        # Pointwise comparison to the analytical Gaussian.
        s_ana = sqrt(var_ana)
        ρ_ana =
            (1 / (s_ana * sqrt(2π))) .*
            exp.(-(grid.centers[1] .- μ_ana) .^ 2 ./ (2 * var_ana))
        L1_err = sum(abs.(ρ_t .- ρ_ana)) * grid.h[1]
        @test L1_err < 5e-3
    end
end

# The Fokker-Planck operator `Qᵀ` and the generator `Q` are matrix
# transposes, so they share the same spectrum (real / complex eigenvalues
# identical, eigenvectors swap left ↔ right).
@testset "Analytical: rate_matrix and fokker_planck_operator share spectrum" begin
    using LinearAlgebra: eigvals
    sys = CoupledSDEs((u, p, t) -> [u[1] - u[1]^3], [0.0]; noise_strength=0.6)
    grid = CartesianGrid((-2.0, 2.0, 100))
    gen = DiffusionGenerator(sys, grid)

    Q = rate_matrix(gen)
    F = fokker_planck_operator(gen)

    λ_Q = sort(real.(eigvals(Matrix(Q))); rev=true)
    λ_F = sort(real.(eigvals(Matrix(F))); rev=true)
    @test maximum(abs.(λ_Q .- λ_F)) < 1e-9
end

# Stationary density of 1D OU `dx = -x dt + σ dW` has variance σ²/2 in
# closed form (Gaussian with that variance is the unique stationary
# solution of the FP equation). Verify quantitatively.
@testset "Analytical: variance of stationary 1D OU = σ²/2" begin
    σ = 1.0
    sys = CoupledSDEs((u, p, t) -> [-u[1]], [0.0]; noise_strength=σ)
    grid = CartesianGrid((-6.0, 6.0, 200))
    gen = DiffusionGenerator(sys, grid)
    ρ = stationary_distribution(gen)
    μ = sum(grid.centers[1] .* ρ) * grid.h[1]
    var_num = sum((grid.centers[1] .- μ) .^ 2 .* ρ) * grid.h[1]
    @test isapprox(μ, 0.0; atol=1e-10)              # symmetric → mean = 0
    @test isapprox(var_num, σ^2 / 2; atol=1e-3)     # var = σ²/2
end

# Propagator semigroup property: exp(t1 · F) · exp(t2 · F) = exp((t1+t2) · F)
# applied to any density. A non-trivial structural test of `propagate_density`.
@testset "Analytical: propagator semigroup property" begin
    sys = CoupledSDEs((u, p, t) -> [-u[1]], [0.0]; noise_strength=1.0)
    grid = CartesianGrid((-6.0, 6.0, 200))
    gen = DiffusionGenerator(sys, grid)

    ρ_0 = zeros(200); ρ_0[100] = 1 / grid.h[1]

    ρs1, _ = propagate_density(gen, 1.5, ρ_0; tol=1e-12)
    ρs2, _ = propagate_density(gen, 1.0, ρs1[:, end]; tol=1e-12)
    ρs_direct, _ = propagate_density(gen, 2.5, ρ_0; tol=1e-12)

    @test maximum(abs.(ρs2[:, end] .- ρs_direct[:, end])) < 1e-10
end

# Eigenmodes of 2D OU `dx = -x dt + σ dW` are -(n+m) for n, m ≥ 0, since
# the operator is separable into two 1D OUs each with Hermite spectrum
# -k for k = 0, 1, 2, .... Multiplicity of eigenvalue -k in 2D is k+1
# (the number of integer pairs (n, m) with n + m = k).
# Slowest 6 (sorted descending): 0, -1 (×2), -2 (×3).
@testset "Analytical: 2D OU eigenmodes are separable sums" begin
    sys = CoupledSDEs((u, p, t) -> [-u[1], -u[2]], [0.0, 0.0]; noise_strength=1.0)
    grid = CartesianGrid((-5.0, 5.0, 71), (-5.0, 5.0, 71))
    gen = DiffusionGenerator(sys, grid)

    λ, _ = eigenmodes(gen, 6)
    expected = [0.0, -1.0, -1.0, -2.0, -2.0, -2.0]
    for i in 1:6
        @test isapprox(real(λ[i]), expected[i]; atol=0.05)
    end
    # Multiplicity-2 pair at -1 is degenerate.
    @test isapprox(real(λ[2]), real(λ[3]); atol=1e-6)
    # Multiplicity-3 cluster at -2.
    @test isapprox(real(λ[4]), real(λ[5]); atol=1e-6)
end

# Reactive rate for pure 1D Brownian motion `dx = σ dW` on (-L, L) with
# reflecting BC and A, B at the leftmost / rightmost cells. With ρ uniform
# = 1/(2L) and linear committor q⁺(x) = (x − x_A)/(x_B − x_A), the rate
# follows from the reversible identity k_AB = ⟨ρ · (σ²/2) · |∇q⁺|²⟩:
#   k_AB = σ² / (2 · 2L · (x_B − x_A)) = σ² / (4 L (x_B − x_A))
@testset "Analytical: reactive rate for pure 1D BM" begin
    σ, L = 1.0, 1.0
    sys = CoupledSDEs((u, p, t) -> [0.0], [0.0]; noise_strength=σ)
    grid = CartesianGrid((-L, L, 401))
    gen = DiffusionGenerator(sys, grid)
    x_A = grid.centers[1][1]
    x_B = grid.centers[1][end]
    A = x -> x[1] < x_A + 1e-12
    B = x -> x[1] > x_B - 1e-12

    res = ReactiveTransition(gen, A, B)
    k_analytic = σ^2 / (2 * 2L * (x_B - x_A))
    @test reactive_rate(res) ≈ k_analytic atol = 1e-5

    # Probability reactive: ∫ ρ q⁺(1−q⁺) dx = (1/(2L)) · (x_B − x_A)/6
    # since ∫₀¹ u(1−u) du = 1/6 after substitution u = (x − x_A)/(x_B − x_A).
    p_analytic = (x_B - x_A) / (12 * L)
    @test probability_reactive(res) ≈ p_analytic atol = 1e-3
end

# =====================================================================
# Type stability
# =====================================================================

# =====================================================================
# Argument validation and edge cases
# =====================================================================

@testset "CartesianGrid: construction validation" begin
    @test_throws ArgumentError CartesianGrid()                            # zero axes
    @test_throws ArgumentError CartesianGrid((1.0, -1.0, 10))             # lo ≥ hi
    @test_throws ArgumentError CartesianGrid((-1.0, 1.0, 1))              # N < 2
    @test_throws ArgumentError CartesianGrid((-1.0, 1.0, 5), (0.0, 0.0, 5))  # second axis lo ≥ hi
end

@testset "DiffusionGenerator: non-diagonal noise rejected" begin
    # Build a CoupledSDEs with a rotated covariance — should be rejected.
    sys = CoupledSDEs((u, p, t) -> [-u[1], -u[2]], [0.0, 0.0]; noise_strength=1.0,
                      covariance=[1.0 0.5; 0.5 1.0])
    grid = CartesianGrid((-2.0, 2.0, 21), (-2.0, 2.0, 21))
    @test_throws ArgumentError DiffusionGenerator(sys, grid)
end

@testset "DiffusionGenerator: dimension mismatch between sys and grid" begin
    sys_1d = CoupledSDEs((u, p, t) -> [-u[1]], [0.0]; noise_strength=1.0)
    grid_2d = CartesianGrid((-2.0, 2.0, 21), (-2.0, 2.0, 21))
    @test_throws DimensionMismatch DiffusionGenerator(sys_1d, grid_2d)
end

@testset "_to_mask: alternative input forms" begin
    sys = CoupledSDEs((u, p, t) -> [-u[1]], [0.0]; noise_strength=1.0)
    grid = CartesianGrid((-2.0, 2.0, 50))
    gen = DiffusionGenerator(sys, grid)

    # Reference predicate
    A_pred = x -> x[1] < -1.0
    qp_pred = forward_committor(gen, A_pred, x -> x[1] > 1.0)

    # BitVector form
    A_bits = BitVector([x < -1.0 for x in grid.centers[1]])
    B_bits = BitVector([x >  1.0 for x in grid.centers[1]])
    qp_bits = forward_committor(gen, A_bits, B_bits)
    @test qp_bits == qp_pred

    # Vector{Int} (linear cell indices)
    A_idx = findall(x -> x < -1.0, grid.centers[1])
    B_idx = findall(x -> x >  1.0, grid.centers[1])
    qp_idx = forward_committor(gen, A_idx, B_idx)
    @test qp_idx == qp_pred

    # Validation: wrong-length BitVector
    @test_throws DimensionMismatch forward_committor(gen, BitVector(falses(20)), B_bits)
    # Validation: out-of-range linear index
    @test_throws BoundsError forward_committor(gen, [60], B_idx)
end

@testset "Committor / MFPT: empty / overlapping sets" begin
    sys = CoupledSDEs((u, p, t) -> [-u[1]], [0.0]; noise_strength=1.0)
    grid = CartesianGrid((-2.0, 2.0, 50))
    gen = DiffusionGenerator(sys, grid)

    @test_throws ArgumentError forward_committor(gen, x -> false, x -> x[1] > 1.0)
    @test_throws ArgumentError forward_committor(gen, x -> x[1] < -1.0, x -> false)
    @test_throws ArgumentError forward_committor(gen, x -> abs(x[1]) < 1.0, x -> abs(x[1]) < 1.0)
    @test_throws ArgumentError backward_committor(gen, x -> false, x -> x[1] > 1.0)
    @test_throws ArgumentError mean_first_passage_time(gen, x -> false)
    @test_throws ArgumentError first_passage_variance(gen, x -> false)
end

@testset "eigenmodes: k clamping and validation" begin
    sys = CoupledSDEs((u, p, t) -> [-u[1]], [0.0]; noise_strength=1.0)
    grid = CartesianGrid((-3.0, 3.0, 30))
    gen = DiffusionGenerator(sys, grid)

    # k = 1 returns just the trivial mode.
    λ1, V1 = eigenmodes(gen, 1)
    @test length(λ1) == 1
    @test size(V1) == (30, 1)
    @test abs(λ1[1]) < 1e-10

    # k > N is clamped to N.
    λ_big, V_big = eigenmodes(gen, 1000)
    @test length(λ_big) == 30
    @test size(V_big) == (30, 30)

    # k < 1 errors.
    @test_throws ArgumentError eigenmodes(gen, 0)
    @test_throws ArgumentError eigenmodes(gen, -1)
end

@testset "CartesianGrid: anisotropic per-axis spacing" begin
    # Stationary density of 2D OU with very different x/y spacings should
    # still match the analytical Gaussian (after normalisation).
    sys = CoupledSDEs((u, p, t) -> [-u[1], -u[2]], [0.0, 0.0]; noise_strength=1.0)
    grid = CartesianGrid((-4.0, 4.0, 81), (-2.0, 2.0, 21))   # h_x ≠ h_y
    gen = DiffusionGenerator(sys, grid)
    @test grid.h[1] != grid.h[2]
    ρ = stationary_distribution(gen)
    xs = [grid.centers[1][I[1]] for I in CartesianIndices(grid.nbox)]
    ys = [grid.centers[2][I[2]] for I in CartesianIndices(grid.nbox)]
    ρ_ana = vec(exp.(-(xs .^ 2 .+ ys .^ 2)) ./ pi)
    @test sum(abs.(ρ .- ρ_ana)) * prod(grid.h) < 1e-2
end

@testset "propagate_density: mass preservation by BC" begin
    sys = CoupledSDEs((u, p, t) -> [-u[1]], [0.0]; noise_strength=1.0)
    grid = CartesianGrid((-3.0, 3.0, 60))

    # Reflecting and Periodic both conserve mass.
    for bc in (Reflecting(), Periodic())
        gen = DiffusionGenerator(sys, grid; bc=bc)
        ρ_0 = zeros(60); ρ_0[30] = 1 / grid.h[1]
        ρs, _ = propagate_density(gen, 5.0, ρ_0; tol=1e-12)
        @test sum(ρs[:, end]) * grid.h[1] ≈ 1.0 atol = 1e-8
    end

    # Absorbing: mass should monotonically decay.
    gen_abs = DiffusionGenerator(sys, grid; bc=Absorbing())
    ρ_0 = zeros(60); ρ_0[30] = 1 / grid.h[1]
    ts = collect(0.0:0.5:5.0)
    ρs, _ = propagate_density(gen_abs, ts[end], ρ_0; Δt=ts[2]-ts[1], tol=1e-12)
    masses = [sum(ρs[:, i]) * grid.h[1] for i in 1:size(ρs, 2)]
    @test issorted(masses; rev=true)            # monotonically decreasing
    @test masses[1] ≈ 1.0 atol = 1e-8           # initial
    @test masses[end] < masses[1]                # actually decayed
end

@testset "Mixed 3D BCs: per-axis tuple in 3D" begin
    sys = CoupledSDEs(
        (u, p, t) -> [-u[1], 0.0, -u[3]], [0.0, 0.0, 0.0]; noise_strength=1.0
    )
    grid = CartesianGrid((-2.0, 2.0, 11), (-pi, pi, 12), (-2.0, 2.0, 11))
    gen = DiffusionGenerator(sys, grid; bc=(Reflecting(), Periodic(), Absorbing()))
    @test gen.bc == (Reflecting(), Periodic(), Absorbing())
    @test gen isa DiffusionGenerator{3,Tuple{Reflecting,Periodic,Absorbing}}

    # Periodic axis preserves row sum on cells where reflecting/periodic are
    # the only contributions; absorbing axis injects negative diagonal on the
    # axis-3 boundary cells. Total row sum is zero only on cells away from
    # absorbing boundaries.
    rsums = vec(sum(gen.Q; dims=2))
    inner_mask = vec([
        (1 < I[3] < grid.nbox[3]) for I in CartesianIndices(grid.nbox)
    ])
    @test maximum(abs, rsums[inner_mask]) < 1e-10
    @test minimum(rsums[.!inner_mask]) < 0       # absorbing axis-3 boundary leaks
end

# =====================================================================
# Grid helpers (non-exported)
# =====================================================================

@testset "grid helpers: not exported but available via :-import" begin
    @test !isdefined(Main, :ball)
    @test !isdefined(Main, :cuboid)
    @test !isdefined(Main, :sublevel)
    @test !isdefined(Main, :reshape_to_grid)
end

@testset "ball / cuboid / sublevel predicates" begin
    using CriticalTransitions: ball, cuboid, sublevel

    A = ball((-1.0, 0.0), 0.25)
    @test A((-1.0, 0.0)) == true
    @test A((-1.2, 0.0)) == true
    @test A((-1.3, 0.0)) == false
    @test A((1.0, 0.0)) == false

    C = cuboid((-1.0, -0.5), (1.0, 0.5))
    @test C((0.0, 0.0)) == true
    @test C((-1.0, -0.5)) == true              # closed box
    @test C((1.0, 0.5)) == true
    @test C((1.1, 0.0)) == false
    @test C((0.0, -0.6)) == false

    f(x) = x[1]^2 + x[2]^2
    S = sublevel(f, 1.0)
    @test S((0.0, 0.0)) == true
    @test S((0.5, 0.5)) == true
    @test S((1.0, 1.0)) == false               # strict <

    # Compose with the public API as a predicate.
    sys = CoupledSDEs((u, p, t) -> [-u[1], -u[2]], [0.0, 0.0]; noise_strength=1.0)
    grid = CartesianGrid((-3.0, 3.0, 21), (-3.0, 3.0, 21))
    gen = DiffusionGenerator(sys, grid)
    qp_ball = forward_committor(gen, A, ball((1.0, 0.0), 0.25))
    @test extrema(qp_ball) == (0.0, 1.0)
end

@testset "reshape_to_grid" begin
    using CriticalTransitions: reshape_to_grid
    sys = CoupledSDEs((u, p, t) -> [-u[1], -u[2]], [0.0, 0.0]; noise_strength=1.0)
    grid = CartesianGrid((-3.0, 3.0, 21), (-3.0, 3.0, 21))
    gen = DiffusionGenerator(sys, grid)

    ρ = stationary_distribution(gen)
    M_gen = reshape_to_grid(ρ, gen)
    M_grid = reshape_to_grid(ρ, grid)
    @test size(M_gen) == (21, 21)
    @test M_gen == M_grid
    @test M_gen[1, 1] == ρ[1]                           # column-major layout

    @test_throws DimensionMismatch reshape_to_grid(zeros(10), grid)
end

@testset "Type stability" begin
    sys = CoupledSDEs((u, p, t) -> [-u[1]], [0.0]; noise_strength=1.0)
    grid = CartesianGrid((-4.0, 4.0, 50))
    gen = DiffusionGenerator(sys, grid)
    @test (@inferred stationary_distribution(gen)) isa Vector{Float64}
    @test (@inferred mean_first_passage_time(gen, x -> abs(x[1]) > 2)) isa Vector{Float64}
end
