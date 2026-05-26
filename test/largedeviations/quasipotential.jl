using CriticalTransitions
using CriticalTransitions: cell_center
using StaticArrays
using LinearAlgebra
using Test

const CT_ = CriticalTransitions

@testset "OLIM quasipotential" begin
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
        @test isapprox(CT_._line_integral(L, SVector(1.0, 0.0), SVector(0.1, 0.0)),
                       0.21; atol = 1.0e-14)

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
        g1 = CartesianGrid((-1.0, 1.0, 10), (-1.0, 1.0, 10))
        g2 = CartesianGrid((-1.0, 1.0, 200), (-1.0, 1.0, 200))
        @test CT_.default_K(g1) == clamp(round(Int, sqrt(10)), 5, 32)
        @test CT_.default_K(g2) == clamp(round(Int, sqrt(200)), 5, 32)
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
            @test isapprox(qp.U[I], dot(x, x); rtol = 0.10)
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

    @testset "back-pointer walk to source" begin
        f(x, p, t) = SVector(-x[1], -x[2])
        sys = CoupledSDEs(f, [0.0, 0.0]; noise_strength = 1.0)
        grid = CartesianGrid((-1.0, 1.0, 31), (-1.0, 1.0, 31))
        qp = quasipotential(sys, grid, [0.0, 0.0]; show_progress = false)
        cur = CartesianIndex(28, 16); visited = [cur]
        while cur != qp.source
            cur = qp.back_pointer[cur].v0
            push!(visited, cur)
            length(visited) > 200 && error("back-pointer walk did not terminate")
        end
        @test cur == qp.source
        # U strictly decreases along the back-pointer chain (Dijkstra invariant).
        Us = [qp.U[I] for I in visited]
        @test all(diff(Us) .<= 1.0e-10)
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
end
