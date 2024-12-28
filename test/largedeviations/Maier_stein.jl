using Test

@testset "Large_deviation Meier Stein" begin
    function meier_stein(u, p, t) # out-of-place
        x, y = u
        dx = x - x^3 - 10 * x * y^2
        dy = -(1 + x^2) * y
        return [dx, dy]
    end
    σ = 0.25
    sys = CoupledSDEs(meier_stein, zeros(2); noise_strength=σ)

    # initial path: parabola
    xx = range(-1.0, 1.0; length=30)
    yy = 0.3 .* (-xx .^ 2 .+ 1)
    init = Matrix([xx yy]')

    x_i = init[:, 1]
    x_f = init[:, end]

    @testset "String method" begin
        string = string_method(
            sys, init; iterations=10_000, ϵ=0.5, show_progress=false
        )
        @test string[:,1] == x_i
        @test string[:,end] == x_f
        @test sum(string[2, :]) ≈ 0 atol=1e-6
        @test sum(string[1, :]) ≈ 0 atol=1e-6
        @test sum(diff(string[1, :])) ≈ 2 atol=1e-6
    end

    @testset "Adam" begin
        gm = geometric_min_action_method(
            sys, x_i, x_f; maxiter=10, verbose=false, show_progress=false
        )
        gm = geometric_min_action_method(
            sys, init; maxiter=500, verbose=false, show_progress=false
        )

        path = gm.path
        action_val = gm.action
        @test all(isapprox.(path[2, :][(end - 5):end], 0, atol=0.01))
        @test all(isapprox.(action_val, 0.3375, atol=0.01))
    end

    @testset "Heteroclinic orbit vs MLP" begin
        import CriticalTransitions as CT
        S(x) = geometric_action(sys, CT.fix_ends(x, init[:, 1], init[:, end]), 1.0)

        gm = geometric_min_action_method(
            sys, init; maxiter=500, verbose=false, show_progress=false
        )
        string = string_method(
            sys, init; iterations=10_000, ϵ=0.5, show_progress=false
        )
        @test S(string) > S(gm.path)
    end
end # gMAM Meier Stein
