using CriticalTransitions
using Test

using CriticalTransitions.CTLibrary: fitzhugh_nagumo

@testset "gMAM FitzHugh-Nagumo" begin
    p = [0.1, 3, 1, 1, 1, 0]
    σ = 0.1
    fhn = CoupledSDEs(fitzhugh_nagumo, zeros(2), p; noise_strength=σ)
    x_i = SA[sqrt(2 / 3), sqrt(2 / 27)]
    x_f = SA[0.001, 0.0]
    res = geometric_min_action_method(fhn, x_i, x_f; N=30, maxiter=500, show_progress=false)
    S = geometric_action(fhn, res.path)
    @test isapprox(S, 0.18, atol=0.01)
end

@testset "gMAM Meier Stein" begin
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

    @testset "LBFGS" begin
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

    @testset "HeymannVandenEijnden" begin # broken
        method = "HeymannVandenEijnden"
        gm = geometric_min_action_method(
            sys, init; maxiter=500, method=method, verbose=false, show_progress=false
        )

        path = gm.path
        action_val = gm.action
        @test all(isapprox.(path[2, :][(end - 5):end], 0, atol=1e-3)) broken = true
        @test all(isapprox.(action_val, 0.3375, atol=1e-3)) broken = true
    end # HeymannVandenEijnden
end # gMAM Meier Stein
