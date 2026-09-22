using CriticalTransitions
using LinearAlgebra
using StaticArrays
using Test

_qp_centers(lo, hi, n) = range(lo + (hi - lo) / (2n), hi - (hi - lo) / (2n); length = n)

@testset "Multiplicative noise end-to-end" begin
    b(u, p, t) = SA[-u[1], -u[2]]
    g(u, p, t) = @SMatrix [1.0 0.0; 0.0 1.0]
    sys = CoupledSDEs(
        b, SA[0.0, 0.0]; g = g,
        noise_prototype = SMatrix{2, 2}(zeros(2, 2)),
    )
    grid = CartesianGrid((-1.0, 1.0, 31), (-1.0, 1.0, 31))
    qp = quasipotential(sys, grid, [0.0, 0.0]; show_progress = false)
    xs = _qp_centers(-1.0, 1.0, 31)
    ys = _qp_centers(-1.0, 1.0, 31)

    @test qp.U[16, 16] == 0.0
    for I in (CartesianIndex(20, 16), CartesianIndex(24, 24), CartesianIndex(12, 20))
        x = SA[xs[I[1]], ys[I[2]]]
        @test isapprox(qp.U[I], dot(x, x); rtol = 0.1)
    end
end

@testset "Regularized quasipotential: equilibrium Langevin" begin
    Vp(x) = x^3 - x
    drift(u, p, t) = SVector(u[2], -Vp(u[1]) - u[2])
    gmat(u, p, t) = @SMatrix [0.0 0.0; 0.0 sqrt(2.0)]
    sys = CoupledSDEs(
        drift, SA[-1.0, 0.0]; g = gmat,
        noise_prototype = SMatrix{2, 2}(zeros(2, 2)),
    )

    grid = CartesianGrid((-1.8, 0.6, 61), (-1.2, 1.2, 61))
    qp = quasipotential(sys, grid, [-1.0, 0.0]; show_progress = false)
    xs = _qp_centers(-1.8, 0.6, 61)
    ps = _qp_centers(-1.2, 1.2, 61)
    isrc = argmin(abs.(xs .+ 1.0))
    jsrc = argmin(abs.(ps))
    isad = argmin(abs.(xs))
    jsad = argmin(abs.(ps))

    @test qp.U[isrc, jsrc] == 0.0
    @test abs(qp.U[isad, jsad] - 0.25) <= 0.02
    @test all(>=(-1.0e-8), filter(isfinite, qp.U))
end

@testset "Regularized quasipotential: non-equilibrium van der Pol" begin
    Vp(x) = x^3 - x
    Dfric(x) = 1.0 - 0.3 * (1 - x^2)
    drift(u, p, t) = SVector(u[2], -Dfric(u[1]) * u[2] - Vp(u[1]))
    gmat(u, p, t) = @SMatrix [0.0 0.0; 0.0 sqrt(2.0)]
    sys = CoupledSDEs(
        drift, SA[-1.0, 0.0]; g = gmat,
        noise_prototype = SMatrix{2, 2}(zeros(2, 2)),
    )
    grid = CartesianGrid((-1.8, 0.6, 41), (-1.2, 1.2, 41))
    qp = quasipotential(sys, grid, [-1.0, 0.0]; show_progress = false)
    xs = _qp_centers(-1.8, 0.6, 41)
    ps = _qp_centers(-1.2, 1.2, 41)
    isrc = argmin(abs.(xs .+ 1.0))
    jsrc = argmin(abs.(ps))
    isad = argmin(abs.(xs))
    jsad = argmin(abs.(ps))

    @test all(isfinite, qp.U)
    @test qp.U[isrc, jsrc] == 0.0
    @test 0.15 < qp.U[isad, jsad] < 0.24
    col = [qp.U[isad, argmin(abs.(ps .- pp))] for pp in (0.0, 0.3, 0.6, 0.9)]
    @test issorted(col)
end
