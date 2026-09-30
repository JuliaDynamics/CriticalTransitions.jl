using Test
using LinearAlgebra: dot
import OptimizationOptimJL

@testset "Ritz method" begin
    @testset "Chebyshev spectral primitives" begin
        n = 8
        nodes = CT._ritz_chebyshev_nodes(n)
        bary = CT._ritz_barycentric_weights(n)
        D = CT._ritz_differentiation_matrix(nodes, bary)

        @test first(nodes) == -1.0
        @test last(nodes) == 1.0
        @test issorted(nodes)

        p = nodes .^ 4 .- 2 .* nodes .^ 2 .+ nodes
        dp = 4 .* nodes .^ 3 .- 4 .* nodes .+ 1
        @test D * p ≈ dp atol = 2.0e-12

        qnodes, qweights = CT._ritz_clenshaw_curtis(40)
        @test sum(qweights) ≈ 2.0 atol = 2.0e-14
        @test dot(qweights, qnodes .^ 4) ≈ 2 / 5 atol = 2.0e-14

        B = CT._ritz_barycentric_matrix(nodes, bary, qnodes)
        C = B * D
        @test B * p ≈ qnodes .^ 4 .- 2 .* qnodes .^ 2 .+ qnodes atol = 2.0e-12
        @test C * p ≈ 4 .* qnodes .^ 3 .- 4 .* qnodes .+ 1 atol = 2.0e-11
    end

    @testset "Known one-dimensional FW action" begin
        drift_1d(u, p, t) = SA[u[1] - u[1]^3]
        sys = CoupledSDEs(drift_1d, SA[-1.0]; noise_strength = 1.0)

        degree = 8
        nodes, _, B, C, qweights = CT._ritz_spectral_matrices(degree, 80, Float64)
        path = reshape((nodes .- 1) ./ 2, 1, :)
        A_at = CT._action_metric(sys)
        S = CT._ritz_on_shell_action(sys, path, B, C, qweights, 0.0, A_at)
        S_E = CT._ritz_on_shell_action(sys, path, B, C, qweights, 0.1, A_at)

        # For b(x) = x - x^3 and the monotone path -1 -> 0,
        # S_0 = -2∫_{-1}^0 b(x) dx = 1/2.
        @test S ≈ 0.5 atol = 2.0e-10
        @test isfinite(S_E)
        @test S_E > S
    end

    @testset "Signed FW energy shell" begin
        constant_drift(u, p, t) = SA[2.0]
        sys = CoupledSDEs(constant_drift, SA[0.0]; noise_strength = 1.0)

        degree = 4
        nodes, _, B, C, qweights = CT._ritz_spectral_matrices(degree, 40, Float64)
        path = reshape((nodes .+ 1) ./ 2, 1, :)
        A_at = CT._action_metric(sys)
        S = CT._ritz_on_shell_action(sys, path, B, C, qweights, -1.0, A_at)

        # Here 2E + |b|^2 = 2 > 0 despite E < 0. The on-shell integrand is constant.
        @test S ≈ 3 / sqrt(2.0) - 2 atol = 2.0e-13
    end

    @testset "Non-positive Onsager-Machlup Lagrangian" begin
        # b(x) = -x, σ = 1 gives
        # L_OM = 1/2 (ẋ + x)^2 - 1/2,
        # which is not positive definite. On the E = 1 shell,
        # |ẋ|^2 = 1 + x^2 and the on-shell action from x=0 to x=1 is
        # 1/2 * (sqrt(2) - asinh(1) + 1).
        linear_drift(u, p, t) = SA[-u[1]]
        sys = CoupledSDEs(linear_drift, SA[0.0]; noise_strength = 1.0)

        degree = 8
        nodes, _, B, C, qweights = CT._ritz_spectral_matrices(degree, 80, Float64)
        path = reshape((nodes .+ 1) ./ 2, 1, :)
        A_at = CT._action_metric(sys)
        S = CT._ritz_on_shell_action(
            sys, path, B, C, qweights, 1.0, A_at, Val(:OM), 1.0
        )
        expected = (sqrt(2.0) - asinh(1.0) + 1.0) / 2

        @test S ≈ expected atol = 2.0e-11
        @test_throws DomainError CT._ritz_on_shell_action(
            sys, path, B, C, qweights, 0.0, A_at, Val(:OM), 1.0
        )
    end

    @testset "Public API" begin
        @test_throws ArgumentError Ritz(degree = 1)
        @test_throws ArgumentError Ritz(degree = 4, quadrature_points = 4)
        @test Ritz(energy = -1.0).energy == -1.0
        @test_throws ArgumentError Ritz(energy = Inf)
        @test_throws ArgumentError Ritz(energy = NaN)

        drift_1d(u, p, t) = SA[u[1] - u[1]^3]
        sys = CoupledSDEs(drift_1d, SA[-1.0]; noise_strength = 1.0)
        method = Ritz(degree = 4, quadrature_points = 40)
        result = minimize_geometric_action(
            sys, SA[-1.0], SA[0.0], method;
            maxiters = 2, output_points = 17, show_progress = false,
        )

        @test result isa MinimumActionPath
        @test length(result.path) == 17
        @test result.path[1] == SA[-1.0]
        @test result.path[end] == SA[0.0]
        @test isfinite(result.action)
        @test result.action ≈ 0.5 atol = 5.0e-4

        @test_throws ArgumentError minimize_geometric_action(
            sys, SA[-1.0], SA[0.0], method;
            functional = "bad", maxiters = 1, show_progress = false,
        )
        @test_throws ArgumentError minimize_geometric_action(
            sys, SA[-1.0], SA[0.0], method;
            functional = "OM", maxiters = 1, show_progress = false,
        )

        linear_drift(u, p, t) = SA[-u[1]]
        om_sys = CoupledSDEs(linear_drift, SA[0.0]; noise_strength = 1.0)
        om_method = Ritz(
            optimizer = OptimizationOptimJL.LBFGS(),
            degree = 4,
            quadrature_points = 40,
            energy = 1.0,
        )
        om_result = minimize_geometric_action(
            om_sys, SA[0.0], SA[1.0], om_method;
            functional = "OM", noise_strength = 1.0, maxiters = 10,
            output_points = 17, show_progress = false,
        )
        expected = (sqrt(2.0) - asinh(1.0) + 1.0) / 2
        @test om_result isa MinimumActionPath
        @test om_result.action ≈ expected atol = 1.0e-8
    end

    @testset "Maier-Stein degree-8 benchmark" begin
        function maier_stein(u, p, t)
            x, y = u
            return SA[x - x^3 - 10 * x * y^2, -(1 + x^2) * y]
        end
        sys = CoupledSDEs(maier_stein, zeros(2); noise_strength = 0.25)

        xx = range(-1.0, 1.0; length = 30)
        yy = 0.3 .* (1 .- xx .^ 2)
        init = Matrix([xx yy]')
        method = Ritz(
            optimizer = OptimizationOptimJL.LBFGS(), degree = 8, quadrature_points = 80
        )
        result = minimize_geometric_action(
            sys, init, method;
            maxiters = 500, abstol = 1.0e-10, reltol = 1.0e-10,
            output_points = 100, show_progress = false,
        )

        @test result.path[1] == SA[-1.0, 0.0]
        @test result.path[end] == SA[1.0, 0.0]
        # Existing gMAM certification for this β=10 Maier-Stein problem is S ≈ 0.3375.
        # The paper's representative Ritz discretization uses degree 8 and nq = 10n.
        @test result.action ≈ 0.3375 atol = 5.0e-3
    end
end
