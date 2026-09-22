using CriticalTransitions
using CriticalTransitions: cell_center
using LinearAlgebra
using StaticArrays
using Test

@testset "Multiplicative noise end-to-end (constant Q)" begin
    b(u, p, t) = SA[-u[1], -u[2]]
    g(u, p, t) = @SMatrix [1.0 0.0; 0.0 1.0]
    sys = CoupledSDEs(
        b, SA[0.0, 0.0]; g = g,
        noise_prototype = SMatrix{2, 2}(zeros(2, 2))
    )
    grid = CartesianGrid((-1.0, 1.0, 31), (-1.0, 1.0, 31))
    qp = quasipotential(sys, grid, [0.0, 0.0]; show_progress = false)
    @test qp.U[16, 16] == 0.0
    for I in (CartesianIndex(20, 16), CartesianIndex(24, 24), CartesianIndex(12, 20))
        x = cell_center(grid, I)
        @test isapprox(qp.U[I], dot(x, x); rtol = 0.1)
    end
end

@testset "regularized OLIM: equilibrium Langevin" begin
    Vp(x) = x^3 - x
    Vpot(x) = x^4 / 4 - x^2 / 2
    Ustar(x, p) = p^2 / 2 + Vpot(x) - Vpot(-1.0)
    drift(u, p, t) = SVector(u[2], -Vp(u[1]) - u[2])
    gmat(u, p, t) = @SMatrix [0.0 0.0; 0.0 sqrt(2.0)]
    sys = CoupledSDEs(drift, SA[-1.0, 0.0]; g = gmat, noise_prototype = SMatrix{2, 2}(zeros(2, 2)))

    grid = CartesianGrid((-1.8, 0.6, 61), (-1.2, 1.2, 61))
    K = 8
    qp = quasipotential(sys, grid, [-1.0, 0.0]; band_radius = K, show_progress = false)
    xs = grid.centers[1]; ps = grid.centers[2]; src = qp.source
    se = Float64[]
    for i in 1:61, j in 1:61
        (-1.0 <= xs[i] <= 0.3 && 0.0 < ps[j] <= 1.0) || continue
        (abs(i - src[1]) <= K && abs(j - src[2]) <= K) && continue
        isfinite(qp.U[i, j]) || continue
        push!(se, (qp.U[i, j] - Ustar(xs[i], ps[j]))^2)
    end
    @test sqrt(sum(se) / length(se)) <= 0.05
    isad = argmin(abs.(xs)); jsad = argmin(abs.(ps))
    @test abs(qp.U[isad, jsad] - 0.25) <= 0.02
    @test qp.U[qp.source] == 0.0
    @test all(>=(-1.0e-8), filter(isfinite, qp.U))
end

@testset "regularized OLIM: van der Pol (non-equilibrium)" begin
    Vp(x) = x^3 - x
    Dfric(x) = 1.0 - 0.3 * (1 - x^2)
    drift(u, p, t) = SVector(u[2], -Dfric(u[1]) * u[2] - Vp(u[1]))
    gmat(u, p, t) = @SMatrix [0.0 0.0; 0.0 sqrt(2.0)]
    sys = CoupledSDEs(drift, SA[-1.0, 0.0]; g = gmat, noise_prototype = SMatrix{2, 2}(zeros(2, 2)))
    grid = CartesianGrid((-1.8, 0.6, 41), (-1.2, 1.2, 41))
    qp = quasipotential(sys, grid, [-1.0, 0.0]; show_progress = false)
    xs = grid.centers[1]; ps = grid.centers[2]
    @test all(isfinite, qp.U)
    @test qp.U[qp.source] == 0.0
    isad = argmin(abs.(xs)); jsad = argmin(abs.(ps))
    @test 0.15 < qp.U[isad, jsad] < 0.24
    col = [qp.U[isad, argmin(abs.(ps .- pp))] for pp in (0.0, 0.3, 0.6, 0.9)]
    @test issorted(col)
end
