using CriticalTransitions
using CriticalTransitions: cell_center

using LinearAlgebra
using Statistics: median
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

        # Passing the endpoint values in must give the identical root while calling `f`
        # two fewer times: that is what lets `_edge_minimum` reuse Φ(0) and Φ(1).
        f = x -> (x - 0.37) * (x + 0.4)
        ncall = Ref(0)
        counted = x -> (ncall[] += 1; f(x))
        rp, okp = CT_.itp_root(counted, 0.0, 1.0; fa = f(0.0), fb = f(1.0))
        withkw = ncall[]
        ncall[] = 0
        rn, okn = CT_.itp_root(counted, 0.0, 1.0)
        @test (rp, okp) === (rn, okn)
        @test withkw == ncall[] - 2
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

    @testset "hermite_U slope limiting" begin
        # The limiter must not perturb data that is already monotone-compatible:
        # secant slopes (the NaN fallback) and zero slopes are inside [0, 3Δ].
        @test CT_._limit_slopes(0.0, 1.0, 1.0, 1.0) == (1.0, 1.0)
        @test CT_._limit_slopes(0.0, 1.0, 0.0, 0.0) == (0.0, 0.0)
        @test CT_._limit_slopes(1.0, 0.0, -1.0, -1.0) == (-1.0, -1.0)
        # Slopes that fight the secant are zeroed; oversteep ones are capped at 3Δ.
        @test CT_._limit_slopes(0.0, 1.0, -5.0, 9.0) == (0.0, 3.0)
        @test CT_._limit_slopes(0.0, 0.0, -1.0, 1.0) == (0.0, 0.0)

        # Regression: unlimited, U0 = U1 = 0 with m0 = -a, m1 = a gives
        # H(λ) = -a λ(1 - λ), i.e. -a/4 at the midpoint. Adding a nonnegative line
        # integral to a negative interpolated value is the only route by which the
        # sweep can return U < 0, so the interpolant must stay in [min(U0,U1), max].
        for λ in (0.1, 0.25, 0.5, 0.75, 0.9)
            @test CT_._hermite_U(0.0, 0.0, -1.0, 1.0, λ) == 0.0
            H = CT_._hermite_U(0.3, 0.7, -4.0, 9.0, λ)
            @test 0.3 - 1.0e-14 <= H <= 0.7 + 1.0e-14
        end
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
        Lf = CT_._GeometricLagrangian{2, Float64}(x -> -x, SMatrix{2, 2, Float64}(I), 1.0e-10)
        @test CT_._diffusion_at(Lf, SVector(0.0, 0.0)) ≈ SMatrix{2, 2, Float64}(I)
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
        @test all(>=(0), filter(isfinite, qp.U))   # exactly 0, not -1e-8
    end

    @testset "rank-1 linear oscillator: U >= 0 and exact CARE reference" begin
        # Linear drift plus constant diffusion has a globally quadratic quasipotential
        # U(x) = x' P x, with P the CARE solution for the *regularized* metric that OLIM
        # actually solves with. That makes this the sharpest available end-to-end check.
        #
        # Regression: at reg = 0.005 (metric anisotropy sqrt(2/reg) = 20) the unlimited
        # Hermite edge interpolant undershot badly -- min U = -0.218 with 123 cells below
        # zero, and U up to 78% *below* the exact value. U ≥ 0 holds by definition and
        # OLIM minimises over a restricted set of paths, so its error must be one-sided
        # (an overshoot); a negative cell or an undershoot is a bug, not discretisation.
        gam = 0.5
        drift(u, p, t) = SVector(u[2], -u[1] - gam * u[2])
        gmat(u, p, t) = @SMatrix [0.0 0.0; 0.0 1.0]
        sys = CoupledSDEs(
            drift, SA[0.0, 0.0]; g = gmat,
            noise_prototype = SMatrix{2, 2}(zeros(2, 2)),
        )
        reg = 0.005
        P = CT_._care([0.0 1.0; -1.0 -gam], [reg 0.0; 0.0 2.0])
        grid = CartesianGrid((-1.5, 1.5, 25), (-1.5, 1.5, 25))
        qp = quasipotential(
            sys, grid, [0.0, 0.0];
            band_radius = 12, regularization = reg, show_progress = false,
        )
        @test all(>=(0), filter(isfinite, qp.U))
        @test qp.U[qp.source] == 0.0
        # Away from the boundary, where the optimal path is not truncated by the box.
        rel = Float64[]
        for I in CartesianIndices(qp.U)
            isfinite(qp.U[I]) || continue
            x = cell_center(grid, I)
            (abs(x[1]) <= 0.9 && abs(x[2]) <= 0.9) || continue
            ue = dot(x, P * x)
            ue > 0.05 && push!(rel, (qp.U[I] - ue) / ue)
        end
        @test minimum(rel) > -0.02                 # was -0.78
        @test median(rel) < 0.1                    # was -0.18 (undershoot)
    end

    @testset "sweep is independent of cell visit order (equivariance)" begin
        # b(u) = (p, -x - γp) is odd and the grid is symmetric about the origin with an
        # odd cell count, so the reflection (x, p) -> (-x, -p) maps cell (i, j) onto
        # (n+1-i, n+1-j) and U must be equivariant. It only can be if every candidate
        # gate is a sound bound: a gate that prunes against the *running* best makes the
        # accepted value depend on the order cells come off the heap, which is not
        # symmetric. This caught pruning that left a 7% asymmetry in a fitted transverse
        # curvature (and up to 8% directly in U) where round-off is the only floor.
        gam = 0.5
        drift(u, p, t) = SVector(u[2], -u[1] - gam * u[2])
        gmat(u, p, t) = @SMatrix [0.0 0.0; 0.0 1.0]
        sys = CoupledSDEs(
            drift, SA[0.0, 0.0]; g = gmat,
            noise_prototype = SMatrix{2, 2}(zeros(2, 2)),
        )
        n = 61
        grid = CartesianGrid((-1.5, 1.5, n), (-1.5, 1.5, n))
        qp = quasipotential(
            sys, grid, [0.0, 0.0];
            band_radius = 8, regularization = 0.05, show_progress = false,
        )
        U = qp.U
        scale = maximum(filter(isfinite, U))
        gaps = Float64[]
        for i in 1:n, j in 1:n
            a, b = U[i, j], U[n + 1 - i, n + 1 - j]
            (isfinite(a) && isfinite(b)) || continue
            m = max(abs(a), abs(b))
            m > 0.01 * scale && push!(gaps, abs(a - b) / m)
        end
        @test maximum(gaps) < 1.0e-6
        @test median(gaps) < 1.0e-9

        # A pop's local updates read only finalised cells and write only non-finalised
        # ones, so batching them across threads cannot change any of them. The field
        # must therefore come out *bitwise* identical, not merely close. Anything that
        # starts reading a `_CONSIDERED` value breaks this rather than degrading it.
        qp_ser = quasipotential(
            sys, grid, [0.0, 0.0];
            band_radius = 8, regularization = 0.05, show_progress = false,
            threaded = false,
        )
        @test all(map(===, qp_ser.U, U))
        @test all(map(===, qp_ser.back_pointer, qp.back_pointer))
    end

    @testset "_sweep! drives a caller-seeded state with no extra keywords" begin
        # Code that needs a custom seed (e.g. a Riccati tube around a limit cycle rather
        # than a point attractor) builds `_OLIMState` itself and hands it to `_sweep!`.
        # Keep that callable with just `verbose`/`show_progress`: the performance
        # keywords added later must stay optional here, not just on `quasipotential`.
        drift(u, p, t) = SVector(-u[1], -u[2])
        sys = CoupledSDEs(drift, SA[0.0, 0.0]; noise_strength = 1.0)
        grid = CartesianGrid((-1.0, 1.0, 21), (-1.0, 1.0, 21))
        L = CT_._geometric_lagrangian(sys, Float64)
        state = CT_._OLIMState(grid, Float64, L.Q_inv isa SMatrix)
        src = CT_._source_cell(grid, SA[0.0, 0.0])
        CT_._sweep!(
            state, grid, src, sys, L, Val(5), Val(3);
            verbose = false, show_progress = false,
        )
        @test all(isfinite, state.U)
        @test all(>=(0), state.U)
        @test state.U[src] == 0
    end

    @testset "incremental update matches the full rescan" begin
        # A pop's only new candidates are the simplexes built from the popped cell;
        # every other front simplex was already offered at an earlier pop. So updating
        # from just those, on top of the value the cell already holds, must reproduce
        # the full rescan. It cannot reproduce it bitwise, because the full rescan also
        # re-evaluates *old* edges whose slope estimates have since improved -- but
        # dropping a candidate can only leave `U` higher, so the difference has a sign.
        gam = 0.5
        drift(u, p, t) = SVector(u[2], -u[1] - gam * u[2])
        gmat(u, p, t) = @SMatrix [0.0 0.0; 0.0 1.0]
        sys = CoupledSDEs(
            drift, SA[0.0, 0.0]; g = gmat,
            noise_prototype = SMatrix{2, 2}(zeros(2, 2)),
        )
        grid = CartesianGrid((-1.5, 1.5, 41), (-1.5, 1.5, 41))
        kw = (;
            band_radius = 8, regularization = 0.02,
            show_progress = false, threaded = false,
        )
        inc = quasipotential(sys, grid, [0.0, 0.0]; kw...)
        full = quasipotential(sys, grid, [0.0, 0.0]; kw..., _full_rescan = true)

        @test all(>=(0), filter(isfinite, inc.U))
        @test count(isfinite, inc.U) == count(isfinite, full.U)
        d = [a - b for (a, b) in zip(inc.U, full.U) if isfinite(a) && isfinite(b)]
        scale = maximum(filter(isfinite, full.U))
        # Dropping candidates can only raise U; allow a hair of slack for the cells
        # where the two sweeps pick different (tied) winners and then diverge.
        @test count(>=(0), d) > 0.99 * length(d)
        @test maximum(d) < 0.01 * scale
        @test -minimum(d) < 1.0e-3 * scale
    end

    @testset "prepared slopes match the per-λ Hermite path" begin
        # `_edge_minimum` limits the slopes once per edge instead of at every λ. That is
        # only sound because `_prepare_slopes` is idempotent: a prepared slope is never
        # NaN and already lies in [0, 3Δ], so re-preparing it is a no-op.
        for (U0, U1, m0, m1) in (
                (0.0, 1.0, NaN, NaN), (1.0, 0.3, 5.0, -2.0), (0.5, 0.5, 0.2, -0.2),
                (0.2, 0.9, 0.1, 4.0), (2.0, 1.0, NaN, 0.7),
            )
            s0, s1 = CT_._prepare_slopes(U0, U1, m0, m1)
            @test CT_._prepare_slopes(U0, U1, s0, s1) === (s0, s1)
            for λ in (0.0, 0.1, 0.5, 0.9, 1.0)
                @test CT_._hermite_prepared(U0, U1, s0, s1, λ) ===
                    CT_._hermite_U(U0, U1, m0, m1, λ)
            end
        end
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
