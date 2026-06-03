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
        Lf = CT_._GeometricLagrangian{2, Float64}(x -> -x, SMatrix{2, 2, Float64}(I), 1e-10)
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
        @test all(>=(-1e-8), filter(isfinite, qp.U))
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
end
