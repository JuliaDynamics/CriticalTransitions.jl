using Test

@testset "Large_deviation Meier Stein" begin
    function meier_stein(u, p, t) # out-of-place
        x, y = u
        dx = x - x^3 - 10 * x * y^2
        dy = -(1 + x^2) * y
        return SA[dx, dy]
    end
    σ = 0.25
    sys = CoupledSDEs(meier_stein, zeros(2); noise_strength = σ)

    # initial path: parabola
    xx = range(-1.0, 1.0; length = 30)
    yy = 0.3 .* (-xx .^ 2 .+ 1)
    init = Matrix([xx yy]')

    x_i = init[:, 1]
    x_f = init[:, end]

    @testset "String method" begin
        for x_init in [init, StateSpaceSet(init')]
            string = string_method(
                sys, x_init; maxiters = 10_000, stepsize = 0.5, show_progress = false
            )
            path = Matrix(string.path)
            @test string.path[1] == x_i
            @test string.path[end] == x_f
            @test sum(path[:, 2]) ≈ 0 atol = 1.0e-6
            @test sum(path[:, 1]) ≈ 0 atol = 1.0e-6
            @test sum(diff(path[:, 1])) ≈ 2 atol = 1.0e-6
        end
    end

    @testset "Adam" begin
        gm = minimize_geometric_action(
            sys, x_i, x_f; maxiters = 10, verbose = false, show_progress = false
        )
        gm = minimize_geometric_action(
            sys, init; maxiters = 500, verbose = false, show_progress = false
        )

        path = Matrix(gm.path)'
        action_val = gm.action
        @test all(isapprox.(path[2, :][(end - 5):end], 0, atol = 0.01))
        @test all(isapprox.(action_val, 0.3375, atol = 0.01))
    end

    @testset "Heteroclinic orbit vs MLP" begin
        import CriticalTransitions as CT
        S(x) = geometric_action(sys, CT.fix_ends(x, init[:, 1], init[:, end]), 1.0)

        gm = minimize_geometric_action(
            sys, init; maxiters = 500, verbose = false, show_progress = false
        )
        string = string_method(
            sys, init; maxiters = 10_000, stepsize = 0.5, show_progress = false
        )
        @test S(permutedims(Matrix(string.path))) > S(Matrix(Matrix(gm.path)'))
    end
end # gMAM Meier Stein

# Heymann-Vanden-Eijnden 2008 §4.1: at β = 1 the Maier-Stein drift is the gradient of
# V(u,v) = -u²/2 + u⁴/4 + (1+u²)v²/2, with minima at (±1, 0) of V = -1/4 and a saddle at
# (0, 0) of V = 0. For additive identity noise the geometric Freidlin-Wentzell action
# from (-1, 0) to (+1, 0) is exactly 2·(V_saddle - V_min) = 1/2.
#
# This is also a regression test for issue #282: an earlier `FW_action` had a spurious
# `/2`, which would make `mlp.action` come back as 1/4 instead of 1/2 here. We additionally
# cross-check that the on-shell quadrature ⟨p, ẋ⟩ used by sgMAM agrees in absolute scale
# with the path-based `fw_action` integrator on the same converged curve.
@testset "Maier-Stein β=1 (gradient): analytic action = 1/2" begin
    function meier_stein_grad(u, p, t)
        x, y = u
        return SA[x - x^3 - x * y^2, -(1 + x^2) * y]
    end
    sys = CoupledSDEs(meier_stein_grad, zeros(2); noise_strength = 1.0)
    N = 100
    xx = range(-1.0, 1.0; length = N)
    yy = 0.3 .* (-xx .^ 2 .+ 1)
    init = Matrix([xx yy]')

    # gMAM (GeometricGradient): step-size 1.0 is fine; backtracking handles stability.
    res_g = minimize_geometric_action(
        sys, init, GeometricGradient(; stepsize = 1.0);
        maxiters = 500, show_progress = false,
    )
    @test isapprox(res_g.action, 0.5; atol = 5.0e-3)

    # sgMAM via Hamiltonian on the same system
    sys_h = FreidlinWentzellHamiltonian(sys)
    res_sg = minimize_simple_geometric_action(
        sys_h, init, GeometricGradient(; stepsize = 1.0);
        maxiters = 500, show_progress = false,
    )
    @test isapprox(res_sg.action, 0.5; atol = 5.0e-3)

    # Issue #282: the sgMAM `mlp.action` (computed via FW_action(xdot, p)) must agree
    # in absolute scale with the path-based fw_action integrator on the same curve.
    # A spurious /2 in FW_action would make this ratio be 1/2 instead of 1.
    sg_path = Matrix(Matrix(res_sg.path)')  # D × N
    T_phys = 1.0  # action is parametrization-invariant on the minimizer
    times = range(0.0, T_phys; length = size(sg_path, 2))
    S_path = fw_action(sys, sg_path, times)
    @test isapprox(res_sg.action, S_path; rtol = 1.0e-2)
end

# Heymann-Vanden-Eijnden 2008 eq. (3.8): "the error on the curve can be made O(Δα²) =
# O(1/N²)". On the smooth β=1 gradient case we have the analytic action S = 1/2, so we
# can check directly that doubling N reduces the action error by roughly 4×.
@testset "gMAM second-order convergence in N (β=1 Maier-Stein)" begin
    function meier_stein_grad(u, p, t)
        x, y = u
        return SA[x - x^3 - x * y^2, -(1 + x^2) * y]
    end
    sys = CoupledSDEs(meier_stein_grad, zeros(2); noise_strength = 1.0)
    actions = Float64[]
    for N in (60, 120)
        xx = range(-1.0, 1.0; length = N)
        yy = 0.3 .* (-xx .^ 2 .+ 1)
        init = Matrix([xx yy]')
        res = minimize_geometric_action(
            sys, init, GeometricGradient(; stepsize = 1.0);
            maxiters = 500, show_progress = false,
        )
        push!(actions, res.action)
    end
    err_coarse = abs(actions[1] - 0.5)
    err_fine = abs(actions[2] - 0.5)
    # Expect ~4× reduction; allow some slack since we're in the asymptotic regime.
    @test err_fine < err_coarse / 2
end
