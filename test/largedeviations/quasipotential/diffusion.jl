using CriticalTransitions
using CriticalTransitions: cell_center
using LinearAlgebra
using StaticArrays
using Test

const CT_ = CriticalTransitions

@testset "Multiplicative noise: Lagrangian and line integral" begin
    # `_QInvDynamic` callable path: build a Lagrangian whose Qinv is a state-
    # dependent SMatrix, then check that values agree with the closed form.
    b = x -> SVector(-x[1], -x[2])
    Qinv_fn = x -> SMatrix{2, 2, Float64}((1 + 0.5 * x[1]^2) * I)
    L_mult = CT_._GeometricLagrangian{2, Float64}(b, Qinv_fn, 1.0e-10)
    x = SVector(0.7, -0.3); v = SVector(0.2, 0.1)
    Qinv = Qinv_fn(x); bx = b(x)
    @test isapprox(
        L_mult(x, v),
        sqrt(dot(v, Qinv * v)) * sqrt(dot(bx, Qinv * bx)) - dot(v, Qinv * bx);
        atol = 1.0e-14,
    )

    # Multiplicative `_line_integral` reduces to Simpson on L(y + s v, v).
    y = SVector(0.5, 0.1)
    v2 = SVector(0.05, -0.02)
    Φ = CT_._line_integral(L_mult, y, v2)
    Φ_simpson = (L_mult(y, v2) + 4 * L_mult(y + 0.5 * v2, v2) + L_mult(y + v2, v2)) / 6
    @test isapprox(Φ, Φ_simpson; atol = 1.0e-14)

    # Consistency: a callable Qinv that returns a *constant* SMatrix must give
    # the same line integral as the additive (SMatrix) path with the same value.
    Qinv_const = SMatrix{2, 2, Float64}(I)
    L_add = CT_._GeometricLagrangian{2, Float64}(b, Qinv_const, 1.0e-10)
    L_dyn = CT_._GeometricLagrangian{2, Float64}(b, _ -> Qinv_const, 1.0e-10)
    @test isapprox(
        CT_._line_integral(L_add, y, v2),
        CT_._line_integral(L_dyn, y, v2);
        atol = 1.0e-14,
    )
end

@testset "Multiplicative noise end-to-end (constant Q)" begin
    # With a *callable* Q(x) that is constant, the result must match the
    # additive solver up to (zero) discretisation difference: it forces the
    # `_QInvDynamic` + multiplicative `_line_integral` code path through the
    # whole sweep, then compares to the analytic U(x) = |x|².
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

@testset "degenerate split" begin
    isdeg, _, _ = CT_._degenerate_split([1.0 0.0; 0.0 1.0])
    @test !isdeg
    isdeg, z, R = CT_._degenerate_split([0.0 0.0; 0.0 2.0])   # noise only in coord 2
    @test isdeg && z == 1 && collect(R) == [2]
    isdeg, z, R = CT_._degenerate_split([0.0 0.0 0.0; 0.0 2.0 0.0; 0.0 0.0 3.0])  # |Z|=1 in 3D
    @test isdeg && z == 1 && collect(R) == [2, 3]
    @test_throws ArgumentError CT_._degenerate_split([0.0 0.0 0.0; 0.0 2.0 0.0; 0.0 0.0 0.0])  # |Z|=2
    @test_throws ArgumentError CT_._degenerate_split([1.0 1.0; 1.0 1.0])  # rotated singular
end

@testset "diffusion accessor" begin
    Lf = CT_._GeometricLagrangian{2, Float64}(x -> -x, SMatrix{2, 2, Float64}(I), 1.0e-10)
    @test CT_._diffusion_at(Lf, SVector(0.0, 0.0)) ≈ SMatrix{2, 2, Float64}(I)
end

@testset "default_regularization" begin
    g80 = CartesianGrid((-1.0, 1.0, 80), (-1.0, 1.0, 80))
    @test CT_.default_regularization(g80) ≈ 0.04
    g160 = CartesianGrid((-1.0, 1.0, 160), (-1.0, 1.0, 160))
    @test CT_.default_regularization(g160) < CT_.default_regularization(g80)
end

@testset "geometric_lagrangian router (regularized rank-1)" begin
    sysf = CoupledSDEs((x, p, t) -> SVector(-x[1], -x[2]), [0.0, 0.0]; noise_strength = 1.0)
    @test CT_._geometric_lagrangian(sysf, Float64) isa CT_._GeometricLagrangian

    drift(u, p, t) = SVector(u[2], -u[1] - u[2])
    gmat(u, p, t) = @SMatrix [0.0 0.0; 0.0 sqrt(2.0)]
    sysd = CoupledSDEs(drift, SA[0.0, 0.0]; g = gmat, noise_prototype = SMatrix{2, 2}(zeros(2, 2)))
    # rank-1 regularized -> ordinary GeometricLagrangian with an invertible (PD) metric
    L = CT_._geometric_lagrangian(sysd, Float64; regularization = 0.04)
    @test L isa CT_._GeometricLagrangian
    @test isposdef(inv(L.Q_inv))
    # rank-1 with no regularization is rejected
    @test_throws ArgumentError CT_._geometric_lagrangian(sysd, Float64)

    drift3(u, p, t) = SVector(u[2], u[3], -u[1])
    g3(u, p, t) = @SMatrix [0.0 0.0 0.0; 0.0 0.0 0.0; 0.0 0.0 1.0]
    sys3 = CoupledSDEs(drift3, SA[0.0, 0.0, 0.0]; g = g3, noise_prototype = SMatrix{3, 3}(zeros(3, 3)))
    @test_throws ArgumentError CT_._geometric_lagrangian(sys3, Float64; regularization = 0.04)

    # state-dependent rank-1 noise -> dynamic regularized metric (_QInvRegDynamic),
    # not the constant SMatrix branch. The regularized metric must be PD everywhere.
    gmul(u, p, t) = @SMatrix [0.0 0.0; 0.0 sqrt(2.0) * (1.0 + 0.2 * u[1]^2)]
    sysm = CoupledSDEs(drift, SA[0.0, 0.0]; g = gmul, noise_prototype = SMatrix{2, 2}(zeros(2, 2)))
    Lm = CT_._geometric_lagrangian(sysm, Float64; regularization = 0.04)
    @test !(Lm.Q_inv isa SMatrix)                                   # dynamic metric
    @test all(isposdef(inv(Lm.Q_inv(SVector(x, p)))) for x in (-0.5, 0.0, 0.5), p in (-0.3, 0.3))
    @test_throws ArgumentError CT_._geometric_lagrangian(sysm, Float64)  # needs reg
end

@testset "regularized OLIM: equilibrium Langevin" begin
    Vp(x) = x^3 - x
    Vpot(x) = x^4 / 4 - x^2 / 2
    Ustar(x, p) = p^2 / 2 + Vpot(x) - Vpot(-1.0)
    drift(u, p, t) = SVector(u[2], -Vp(u[1]) - u[2])         # gamma = 1, closed form U*
    gmat(u, p, t) = @SMatrix [0.0 0.0; 0.0 sqrt(2.0)]
    sys = CoupledSDEs(drift, SA[-1.0, 0.0]; g = gmat, noise_prototype = SMatrix{2, 2}(zeros(2, 2)))

    grid = CartesianGrid((-1.8, 0.6, 61), (-1.2, 1.2, 61))
    qp = quasipotential(sys, grid, [-1.0, 0.0]; show_progress = false)
    xs = grid.centers[1]; ps = grid.centers[2]; src = qp.source
    K = CT_.default_K(grid)
    # escape sheet (p > 0): the regularization bias is negligible there
    se = Float64[]
    for i in 1:61, j in 1:61
        (-1.0 <= xs[i] <= 0.3 && 0.0 < ps[j] <= 1.0) || continue
        (abs(i - src[1]) <= K && abs(j - src[2]) <= K) && continue
        isfinite(qp.U[i, j]) || continue
        push!(se, (qp.U[i, j] - Ustar(xs[i], ps[j]))^2)
    end
    @test sqrt(sum(se) / length(se)) <= 0.05                 # escape-sheet RMS
    isad = argmin(abs.(xs)); jsad = argmin(abs.(ps))
    @test abs(qp.U[isad, jsad] - 0.25) <= 0.02               # saddle barrier
    @test qp.U[qp.source] == 0.0
    @test all(>=(-1.0e-8), filter(isfinite, qp.U))
end

@testset "regularized OLIM: van der Pol (non-equilibrium)" begin
    Vp(x) = x^3 - x
    Dfric(x) = 1.0 - 0.3 * (1 - x^2)                         # state-dependent friction
    drift(u, p, t) = SVector(u[2], -Dfric(u[1]) * u[2] - Vp(u[1]))
    gmat(u, p, t) = @SMatrix [0.0 0.0; 0.0 sqrt(2.0)]
    sys = CoupledSDEs(drift, SA[-1.0, 0.0]; g = gmat, noise_prototype = SMatrix{2, 2}(zeros(2, 2)))
    grid = CartesianGrid((-1.8, 0.6, 41), (-1.2, 1.2, 41))
    qp = quasipotential(sys, grid, [-1.0, 0.0]; show_progress = false)
    xs = grid.centers[1]; ps = grid.centers[2]
    @test all(isfinite, qp.U)                                # full field, no dead band
    @test qp.U[qp.source] == 0.0
    isad = argmin(abs.(xs)); jsad = argmin(abs.(ps))
    @test 0.15 < qp.U[isad, jsad] < 0.24                     # finite, below equilibrium 0.25
    col = [qp.U[isad, argmin(abs.(ps .- pp))] for pp in (0.0, 0.3, 0.6, 0.9)]
    @test issorted(col)                                      # monotone escape sheet
end
