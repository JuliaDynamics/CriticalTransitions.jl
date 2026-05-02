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
# Type stability
# =====================================================================

@testset "Type stability" begin
    sys = CoupledSDEs((u, p, t) -> [-u[1]], [0.0]; noise_strength=1.0)
    grid = CartesianGrid((-4.0, 4.0, 50))
    gen = DiffusionGenerator(sys, grid)
    @test (@inferred stationary_distribution(gen)) isa Vector{Float64}
    @test (@inferred mean_first_passage_time(gen, x -> abs(x[1]) > 2)) isa Vector{Float64}
end
