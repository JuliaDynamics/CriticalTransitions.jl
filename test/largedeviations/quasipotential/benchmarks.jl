using CriticalTransitions
using CriticalTransitions: cell_center
using LinearAlgebra
using StaticArrays
using Test

@testset "3D K_seed=0 bootstrap" begin
    # Regression: pre-fix, K_seed=0 in 3D left every non-source cell at Inf
    # because the 3D simplex pass never produced a standalone Φ0 vertex.
    f(x, p, t) = SA[-x[1], -x[2], -x[3]]
    sys = CoupledSDEs(f, [0.0, 0.0, 0.0]; noise_strength = 1.0)
    grid = CartesianGrid((-1.0, 1.0, 9), (-1.0, 1.0, 9), (-1.0, 1.0, 9))
    qp = quasipotential(
        sys, grid, [0.0, 0.0, 0.0];
        show_progress = false, near_source_layers = 0, band_radius = 2,
    )
    @test all(isfinite, qp.U)
    @test qp.U[5, 5, 5] == 0.0
    @test all(>=(0), qp.U)
end

@testset "2D quadratic well end-to-end" begin
    f(x, p, t) = SVector(-x[1], -x[2])
    sys = CoupledSDEs(f, [0.0, 0.0]; noise_strength = 1.0)
    grid = CartesianGrid((-1.0, 1.0, 31), (-1.0, 1.0, 31))
    qp = quasipotential(sys, grid, [0.0, 0.0]; show_progress = false)
    @test qp.U[16, 16] == 0.0  # source cell
    # Analytic: gradient drift b = -∇V with V = |x|²/2; FW quasipotential U = 2 V = |x|².
    # Check at several cells along the diagonal and on-axis.
    for I in (CartesianIndex(20, 16), CartesianIndex(24, 24), CartesianIndex(12, 20))
        x = cell_center(grid, I)
        @test isapprox(qp.U[I], dot(x, x); rtol = 0.1)
    end
    # U(x) ≥ 0 by definition; -ε is a real bug, not just an interpolation artefact.
    @test all(>=(0), filter(isfinite, qp.U))
end

@testset "D=5 warning" begin
    f(x, p, t) = -x
    sys = CoupledSDEs(f, zeros(5); noise_strength = 1.0)
    grid = CartesianGrid(
        (-1.0, 1.0, 5), (-1.0, 1.0, 5), (-1.0, 1.0, 5),
        (-1.0, 1.0, 5), (-1.0, 1.0, 5),
    )
    @test_logs (:warn, r"D=5") quasipotential(
        sys, grid, zeros(5); band_radius = 3,
        near_source_layers = 0, show_progress = false,
    )
end

@testset "3D gradient well end-to-end" begin
    # D=3 exercises `_add_simplex_candidates{3}` (triangle Newton) and
    # `_triangle_minimum`. Analytic: b = -∇V with V = |x|²/2, so U = 2V = |x|².
    f(x, p, t) = SVector(-x[1], -x[2], -x[3])
    sys = CoupledSDEs(f, [0.0, 0.0, 0.0]; noise_strength = 1.0)
    grid = CartesianGrid((-1.0, 1.0, 11), (-1.0, 1.0, 11), (-1.0, 1.0, 11))
    qp = quasipotential(
        sys, grid, [0.0, 0.0, 0.0];
        show_progress = false, band_radius = 3
    )
    @test qp.U[6, 6, 6] == 0.0  # source cell at origin
    for I in (
            CartesianIndex(9, 6, 6),
            CartesianIndex(9, 9, 6),
            CartesianIndex(9, 9, 9),
        )
        x = cell_center(grid, I)
        @test isapprox(qp.U[I], dot(x, x); rtol = 0.1)
    end
    @test all(>=(0), filter(isfinite, qp.U))
end

@testset "Maier-Stein non-gradient" begin
    f(x, p, t) = SVector(
        x[1] - x[1]^3 - 5 * x[1] * x[2]^2,
        -(1 + x[1]^2) * x[2],
    )
    sys = CoupledSDEs(f, [-1.0, 0.0]; noise_strength = 0.3)
    grid = CartesianGrid((-1.5, 1.5, 61), (-1.0, 1.0, 41))
    qp = quasipotential(sys, grid, [-1.0, 0.0]; show_progress = false)
    # Saddle at (0, 0) maps to grid cell (31, 21). FW quasipotential barrier
    # from (-1, 0) to the saddle for this Maier-Stein system is U_saddle ≈ 0.5.
    saddle = CartesianIndex(31, 21)
    @test isapprox(qp.U[saddle], 0.5; rtol = 0.15)
    # The two stable fixed points are symmetric under x → -x; both attractors
    # should be reached by the sweep with U ≥ 0 throughout.
    @test all(>=(0), filter(isfinite, qp.U))
end
