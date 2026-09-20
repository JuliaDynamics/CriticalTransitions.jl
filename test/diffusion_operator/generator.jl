using CriticalTransitions
using DoubleFloats: Double64
using LinearAlgebra: diag
using SparseArrays: rowvals, nonzeros, nzrange, SparseMatrixCSC, sparse
using Test

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

@testset "DiffusionGenerator (1D OU): rate-matrix structure" begin
    sys = CoupledSDEs((u, p, t) -> [-u[1]], [0.0]; noise_strength = 1.0)
    grid = CartesianGrid((-4.0, 4.0, 200))
    gen = DiffusionGenerator(sys, grid)
    @test gen isa DiffusionGenerator
    Q = gen.Q
    @test Q isa SparseMatrixCSC{Float64, Int}
    @test size(Q) == (200, 200)
    @test maximum(abs, vec(sum(Q; dims = 2))) < 1.0e-12
    omin, omax = _offdiag_extrema(Q)
    @test omin >= 0           # rate matrix: off-diagonals = transition rates ≥ 0
    @test isfinite(omax)
    @test all(diag(Q) .<= 0)  # rate matrix: diagonal = -escape rate ≤ 0
end

@testset "rate_matrix and m_matrix accessors" begin
    sys = CoupledSDEs((u, p, t) -> [-u[1]], [0.0]; noise_strength = 1.0)
    grid = CartesianGrid((-4.0, 4.0, 200))
    gen = DiffusionGenerator(sys, grid)
    Q = rate_matrix(gen)
    L = m_matrix(gen)
    @test Q === gen.Q                              # alias for the field
    @test maximum(abs, L + Q) < 1.0e-12              # M-matrix is the negation
    omin, omax = _offdiag_extrema(L)
    @test omax <= 0           # M-matrix: off-diagonals ≤ 0
    @test all(diag(L) .>= 0)  # M-matrix: diagonal ≥ 0
end

@testset "fokker_planck_operator" begin
    sys = CoupledSDEs((u, p, t) -> [-u[1]], [0.0]; noise_strength = 1.0)
    grid = CartesianGrid((-4.0, 4.0, 50))
    gen = DiffusionGenerator(sys, grid)
    F = fokker_planck_operator(gen)
    @test F == sparse(transpose(gen.Q))
    # Stationary density solves F ρ = 0 (FP nullspace).
    ρ = stationary_distribution(gen)
    @test maximum(abs, F * ρ) < 1.0e-10
end

@testset "BC: Reflecting default" begin
    sys = CoupledSDEs((u, p, t) -> [-u[1]], [0.0]; noise_strength = 1.0)
    grid = CartesianGrid((-4.0, 4.0, 50))
    gen = DiffusionGenerator(sys, grid)
    @test gen.bc == (Reflecting(),)
    @test gen isa DiffusionGenerator{1, Tuple{Reflecting}}
    @test maximum(abs, vec(sum(gen.Q; dims = 2))) < 1.0e-12
end

@testset "BC: Periodic — 1D ring with pure diffusion" begin
    sys = CoupledSDEs((u, p, t) -> [0.0], [0.0]; noise_strength = 1.0)
    grid = CartesianGrid((-pi, pi, 60))
    gen = DiffusionGenerator(sys, grid; bc = Periodic())
    @test gen.bc == (Periodic(),)
    @test gen isa DiffusionGenerator{1, Tuple{Periodic}}
    @test maximum(abs, vec(sum(gen.Q; dims = 2))) < 1.0e-12
    ρ = stationary_distribution(gen)
    @test maximum(ρ) - minimum(ρ) < 1.0e-12
end

@testset "BC: Absorbing leaks mass on boundary" begin
    sys = CoupledSDEs((u, p, t) -> [-u[1]], [0.0]; noise_strength = 1.0)
    grid = CartesianGrid((-3.0, 3.0, 60))
    gen = DiffusionGenerator(sys, grid; bc = Absorbing())
    @test gen.bc == (Absorbing(),)
    rsums = vec(sum(gen.Q; dims = 2))
    @test maximum(abs, rsums[2:59]) < 1.0e-10
    @test rsums[1] < -1.0e-3
    @test rsums[end] < -1.0e-3
end

@testset "BC: per-axis tuple" begin
    sys = CoupledSDEs((u, p, t) -> [-u[1], 0.0], [0.0, 0.0]; noise_strength = 1.0)
    grid = CartesianGrid((-2.0, 2.0, 21), (-pi, pi, 24))
    gen = DiffusionGenerator(sys, grid; bc = (Reflecting(), Periodic()))
    @test gen.bc == (Reflecting(), Periodic())
    @test gen isa DiffusionGenerator{2, Tuple{Reflecting, Periodic}}
    @test maximum(abs, vec(sum(gen.Q; dims = 2))) < 1.0e-12
end

@testset "BC: validation errors" begin
    sys = CoupledSDEs((u, p, t) -> [-u[1]], [0.0]; noise_strength = 1.0)
    grid = CartesianGrid((-3.0, 3.0, 20))
    # Symbols no longer accepted (catches typos at the type system).
    @test_throws ArgumentError DiffusionGenerator(sys, grid; bc = :reflecting)
    # Wrong tuple length.
    @test_throws ArgumentError DiffusionGenerator(
        sys, grid; bc = (Reflecting(), Periodic())
    )
end

@testset "BC: stationary_distribution rejects Absorbing" begin
    sys = CoupledSDEs((u, p, t) -> [-u[1]], [0.0]; noise_strength = 1.0)
    grid = CartesianGrid((-3.0, 3.0, 30))
    gen_abs = DiffusionGenerator(sys, grid; bc = Absorbing())
    @test_throws ArgumentError stationary_distribution(gen_abs)
end

@testset "CartesianGrid: construction validation" begin
    @test_throws ArgumentError CartesianGrid()                            # zero axes
    @test_throws ArgumentError CartesianGrid((1.0, -1.0, 10))             # lo ≥ hi
    @test_throws ArgumentError CartesianGrid((-1.0, 1.0, 1))              # N < 2
    @test_throws ArgumentError CartesianGrid((-1.0, 1.0, 5), (0.0, 0.0, 5))  # second axis lo ≥ hi
end

@testset "DiffusionGenerator: non-diagonal noise rejected" begin
    # Build a CoupledSDEs with a rotated covariance — should be rejected.
    sys = CoupledSDEs(
        (u, p, t) -> [-u[1], -u[2]], [0.0, 0.0]; noise_strength = 1.0,
        covariance = [1.0 0.5; 0.5 1.0]
    )
    grid = CartesianGrid((-2.0, 2.0, 21), (-2.0, 2.0, 21))
    @test_throws ArgumentError DiffusionGenerator(sys, grid)
end

@testset "DiffusionGenerator: dimension mismatch between sys and grid" begin
    sys_1d = CoupledSDEs((u, p, t) -> [-u[1]], [0.0]; noise_strength = 1.0)
    grid_2d = CartesianGrid((-2.0, 2.0, 21), (-2.0, 2.0, 21))
    @test_throws DimensionMismatch DiffusionGenerator(sys_1d, grid_2d)
end

@testset "CartesianGrid: anisotropic per-axis spacing" begin
    # Stationary density of 2D OU with very different x/y spacings should
    # still match the analytical Gaussian (after normalisation).
    sys = CoupledSDEs((u, p, t) -> [-u[1], -u[2]], [0.0, 0.0]; noise_strength = 1.0)
    grid = CartesianGrid((-4.0, 4.0, 81), (-2.0, 2.0, 21))   # h_x ≠ h_y
    gen = DiffusionGenerator(sys, grid)
    @test grid.h[1] != grid.h[2]
    ρ = stationary_distribution(gen)
    xs = [grid.centers[1][I[1]] for I in CartesianIndices(grid.nbox)]
    ys = [grid.centers[2][I[2]] for I in CartesianIndices(grid.nbox)]
    ρ_ana = vec(exp.(-(xs .^ 2 .+ ys .^ 2)) ./ pi)
    @test sum(abs.(ρ .- ρ_ana)) * prod(grid.h) < 1.0e-2
end

@testset "Mixed 3D BCs: per-axis tuple in 3D" begin
    sys = CoupledSDEs(
        (u, p, t) -> [-u[1], 0.0, -u[3]], [0.0, 0.0, 0.0]; noise_strength = 1.0
    )
    grid = CartesianGrid((-2.0, 2.0, 11), (-pi, pi, 12), (-2.0, 2.0, 11))
    gen = DiffusionGenerator(sys, grid; bc = (Reflecting(), Periodic(), Absorbing()))
    @test gen.bc == (Reflecting(), Periodic(), Absorbing())
    @test gen isa DiffusionGenerator{3, Tuple{Reflecting, Periodic, Absorbing}}

    # Periodic axis preserves row sum on cells where reflecting/periodic are
    # the only contributions; absorbing axis injects negative diagonal on the
    # axis-3 boundary cells. Total row sum is zero only on cells away from
    # absorbing boundaries.
    rsums = vec(sum(gen.Q; dims = 2))
    inner_mask = vec(
        [
            (1 < I[3] < grid.nbox[3]) for I in CartesianIndices(grid.nbox)
        ]
    )
    @test maximum(abs, rsums[inner_mask]) < 1.0e-10
    @test minimum(rsums[.!inner_mask]) < 0       # absorbing axis-3 boundary leaks
end

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
    sys = CoupledSDEs((u, p, t) -> [-u[1], -u[2]], [0.0, 0.0]; noise_strength = 1.0)
    grid = CartesianGrid((-3.0, 3.0, 21), (-3.0, 3.0, 21))
    gen = DiffusionGenerator(sys, grid)
    target = ball((1.0, 0.0), 0.25)
    τ = mean_first_passage_time(gen, target)
    @test minimum(τ) >= 0
end

@testset "reshape_to_grid" begin
    using CriticalTransitions: reshape_to_grid
    sys = CoupledSDEs((u, p, t) -> [-u[1], -u[2]], [0.0, 0.0]; noise_strength = 1.0)
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

@testset "Parametric float type: CartesianGrid{T}/DiffusionGenerator{T}" begin
    sys = CoupledSDEs((u, p, t) -> [-u[1]], [0.0]; noise_strength = 1.0)

    # Default constructor remains Float64.
    grid_f = CartesianGrid((-4.0, 4.0, 60))
    @test grid_f isa CartesianGrid{1, Float64}
    @test CriticalTransitions.floattype(grid_f) === Float64
    gen_f = DiffusionGenerator(sys, grid_f)
    @test gen_f isa DiffusionGenerator{1, Tuple{Reflecting}, Float64}
    @test eltype(gen_f.Q) === Float64

    # Explicit Double64 grid → Double64 generator, Double64 matrix entries.
    grid_d = CartesianGrid{Double64}((-4.0, 4.0, 60))
    @test grid_d isa CartesianGrid{1, Double64}
    @test CriticalTransitions.floattype(grid_d) === Double64
    @test eltype(grid_d.h) === Double64
    gen_d = DiffusionGenerator(sys, grid_d)
    @test gen_d isa DiffusionGenerator{1, Tuple{Reflecting}, Double64}
    @test eltype(gen_d.Q) === Double64

    # Same algebra at higher precision: Q at Double64 differs from Q at
    # Float64 by a few hundred eps(Float64) (a few ulps per matrix entry,
    # accumulated across the SG stencil). Densify both sides first; sparse
    # broadcasting across mixed eltypes returns SparseMatrixCSC{Any} on
    # Julia 1.10 and trips `zero(Any)` in `maximum`.
    Q_f_d = Matrix{Double64}(gen_f.Q)
    Q_d = Matrix(gen_d.Q)
    @test Float64(maximum(abs.(Q_f_d .- Q_d))) < 1.0e-13

    # Mass conservation holds at Double64 precision.
    @test Float64(maximum(abs, vec(sum(gen_d.Q; dims = 2)))) < 1.0e-25

    # Float64 solve path still works through the generic interface.
    ρ_f = stationary_distribution(gen_f)
    @test ρ_f isa Vector{Float64}
    @test sum(ρ_f) * CriticalTransitions.cell_volume(grid_f) ≈ 1.0 atol = 1.0e-10
end
