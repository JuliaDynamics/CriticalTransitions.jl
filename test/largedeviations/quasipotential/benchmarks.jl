using CriticalTransitions
using LinearAlgebra
using StaticArrays
using Test

_grid_coordinate(lo, hi, n, i) = lo + (i - 0.5) * (hi - lo) / n

@testset "3D K_seed=0 bootstrap" begin
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
    @test qp.U[16, 16] == 0.0
    for I in (CartesianIndex(20, 16), CartesianIndex(24, 24), CartesianIndex(12, 20))
        x = SA[
            _grid_coordinate(-1.0, 1.0, 31, I[1]),
            _grid_coordinate(-1.0, 1.0, 31, I[2]),
        ]
        @test isapprox(qp.U[I], dot(x, x); rtol = 0.1)
    end
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
    f(x, p, t) = SVector(-x[1], -x[2], -x[3])
    sys = CoupledSDEs(f, [0.0, 0.0, 0.0]; noise_strength = 1.0)
    grid = CartesianGrid((-1.0, 1.0, 11), (-1.0, 1.0, 11), (-1.0, 1.0, 11))
    qp = quasipotential(
        sys, grid, [0.0, 0.0, 0.0];
        show_progress = false, band_radius = 3,
    )
    @test qp.U[6, 6, 6] == 0.0
    for I in (
            CartesianIndex(9, 6, 6),
            CartesianIndex(9, 9, 6),
            CartesianIndex(9, 9, 9),
        )
        x = SA[
            _grid_coordinate(-1.0, 1.0, 11, I[1]),
            _grid_coordinate(-1.0, 1.0, 11, I[2]),
            _grid_coordinate(-1.0, 1.0, 11, I[3]),
        ]
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
    saddle = CartesianIndex(31, 21)
    @test isapprox(qp.U[saddle], 0.5; rtol = 0.15)
    @test all(>=(0), filter(isfinite, qp.U))
end
