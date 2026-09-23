using CriticalTransitions
using LinearAlgebra
using StaticArrays
using Test

const CT_ = CriticalTransitions

@testset "itp_root" begin
    r, ok = CT_.itp_root(x -> (x - 0.37) * (x + 0.4), 0.0, 1.0)
    @test ok
    @test isapprox(r, 0.37; atol = 1.0e-10)

    r2, ok2 = CT_.itp_root(x -> x + 1.0, 0.0, 1.0)
    @test !ok2

    rf, _ = CT_.itp_root(x -> x - 0.25f0, 0.0f0, 1.0f0)
    @test isapprox(rf, 0.25f0; atol = 1.0f-5)
end

@testset "GeometricLagrangian" begin
    Qinv = SMatrix{2, 2, Float64}(I)
    L = CT_._GeometricLagrangian{2, Float64}(x -> -x, Qinv, 1.0e-10)
    x = SVector(1.0, 0.5); v = SVector(0.2, -0.1)
    bx = -x
    @test isapprox(
        L(x, v), sqrt(dot(v, v)) * sqrt(dot(bx, bx)) - dot(v, bx);
        atol = 1.0e-14,
    )

    L0 = CT_._GeometricLagrangian{2, Float64}(x -> zero(x), Qinv, 1.0e-10)
    @test L0(SVector(0.0, 0.0), SVector(0.1, 0.2)) ==
        0.5 * dot(SVector(0.1, 0.2), SVector(0.1, 0.2))
end

@testset "line_integral / hermite_U" begin
    Qinv = SMatrix{2, 2, Float64}(I)
    L = CT_._GeometricLagrangian{2, Float64}(x -> -x, Qinv, 1.0e-10)
    # L_g(x, v) = |v||x| + v·x along y + s v with y=(1,0), v=(0.1,0), x(s) = (1+0.1s, 0)
    # gives 0.2*(1 + 0.1 s), exact integral 0.21 (Simpson is exact on linear integrands).
    @test isapprox(
        CT_._line_integral(L, SVector(1.0, 0.0), SVector(0.1, 0.0)),
        0.21; atol = 1.0e-14
    )

    @test isapprox(CT_._hermite_U(0.0, 1.0, 0.0, 0.0, 0.5), 0.5; atol = 1.0e-12)
    @test isapprox(
        CT_._hermite_U(0.0, 1.0, 0.0, 0.0, 0.25),
        3 * 0.25^2 - 2 * 0.25^3;
        atol = 1.0e-12,
    )
    # NaN slopes fall back to the secant (U1 - U0), reducing to linear interpolation.
    @test isapprox(CT_._hermite_U(0.0, 1.0, NaN, NaN, 0.5), 0.5; atol = 1.0e-12)
end

@testset "stencil_offsets" begin
    s2 = CT_._stencil_offsets(Val(3), Val(2))
    @test all(o -> 0 < o[1]^2 + o[2]^2 <= 9, s2)
    @test length(s2) == 28  # integer lattice points 0 < |δ|² ≤ 9 in ℤ²

    s3 = CT_._stencil_offsets(Val(2), Val(3))
    @test all(o -> 0 < o[1]^2 + o[2]^2 + o[3]^2 <= 4, s3)
    @test length(s3) == 32  # integer lattice points 0 < |δ|² ≤ 4 in ℤ³
end

@testset "OLIMState initialisation" begin
    grid = CartesianGrid((-1.0, 1.0, 10), (-1.0, 1.0, 10))
    st = CT_._OLIMState(grid, Float64)
    @test all(isinf, st.U)
    @test all(==(CT_._UNKNOWN), st.status)
    @test count(st.front) == 0
    @test all(==(CT_.BackRef{2}()), st.back_pointer)
end

@testset "default_K" begin
    g2_small = CartesianGrid((-1.0, 1.0, 10), (-1.0, 1.0, 10))
    g2_large = CartesianGrid((-1.0, 1.0, 200), (-1.0, 1.0, 200))
    @test CT_.default_K(g2_small) == 5            # 2D floor
    @test CT_.default_K(g2_large) == 14           # sqrt(200) ≈ 14, < cap 32
    g3_small = CartesianGrid((-1.0, 1.0, 5), (-1.0, 1.0, 5), (-1.0, 1.0, 5))
    g3_mid = CartesianGrid((-1.0, 1.0, 9), (-1.0, 1.0, 9), (-1.0, 1.0, 9))
    g3_large = CartesianGrid((-1.0, 1.0, 200), (-1.0, 1.0, 200), (-1.0, 1.0, 200))
    @test CT_.default_K(g3_small) == 2            # 3D floor
    @test CT_.default_K(g3_mid) == 3            # sqrt(9) = 3
    @test CT_.default_K(g3_large) == 8            # 3D cap
    g4 = CartesianGrid((-1.0, 1.0, 5), (-1.0, 1.0, 5), (-1.0, 1.0, 5), (-1.0, 1.0, 5))
    @test CT_.default_K(g4) == 2                  # 4D+ floor
end

@testset "cached line integral == live (additive)" begin
    # The cached additive path must be bitwise-identical to the live path; this
    # locks the equivalence the drift cache relies on.
    b = x -> SVector(-x[1], x[1] - x[2])
    Qinv = SMatrix{2, 2, Float64}(2.0, 0.3, 0.3, 1.5)   # symmetric PD
    L = CT_._GeometricLagrangian{2, Float64}(b, Qinv, 1.0e-10)
    nd(x) = (bb = b(x); q = dot(bb, Qinv * bb); CT_._NodeData{2, Float64}(bb, q, sqrt(q)))
    y = SVector(0.2, -0.4); v = SVector(0.1, 0.25)
    live = CT_._line_integral(L, y, v)
    @test CT_._line_integral(L, y, v, nd(y), nd(y + v)) == live   # both endpoints cached
    @test CT_._line_integral(L, y, v, nothing, nd(y + v)) == live # only s=1 cached (edge use)
    @test CT_._line_integral(L, y, v, nd(y), nothing) == live     # only s=0 cached
    # near-zero-drift branch (|b|²_Q < eps_b²) must also match
    bz = x -> SVector(0.0, 0.0)
    Lz = CT_._GeometricLagrangian{2, Float64}(bz, Qinv, 1.0e-10)
    ndz = CT_._NodeData{2, Float64}(SVector(0.0, 0.0), 0.0, 0.0)
    @test CT_._line_integral(Lz, y, v, ndz, ndz) == CT_._line_integral(Lz, y, v)
end
