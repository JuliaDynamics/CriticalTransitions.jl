@testset "gMAM FitzHugh-Nagumo" begin
    p = [0.1, 3, 1, 1, 1, 0]; σ = 0.1
    fhn = CoupledSDEs(fitzhugh_nagumo, idfunc, zeros(2), σ, p)
    x_i = SA[sqrt(2 / 3), sqrt(2 / 27)]
    x_f = SA[0.001, 0.0]
    N = 100
    res = geometric_min_action_method(fhn, x_i, x_f; N=75, maxiter=200)
    S = geometric_action(fhn, res[1][end])
    @test isapprox(S, 0.18, atol=0.01)
end

"""
@testset "gMAM Meier Stein" begin
    function meier_stein(u, p, t) # out-of-place
        x, y = u
        dx = x - x^3 - 10 * x * y^2
        dy = -(1 + x^2) * y
        [dx, dy]
    end
    σ = 0.25
    sys = StochSystem(meier_stein, [], zeros(2), σ, idfunc, nothing, I(2), "WhiteGauss")

    # initial path: parabola
    xx = range(-1.0, 1.0, length = 30)
    yy = 0.3 .* (-xx .^ 2 .+ 1)
    init = Matrix([xx yy]')

    x_i = init[:, 1]
    x_f = init[:, end]

    @testset "LBFGS" begin
        gm = geometric_min_action_method(sys, x_i, x_f, maxiter = 10, verbose = false)# runtest
        gm = geometric_min_action_method(sys, init, maxiter = 100, verbose=false)

        path = gm[1][end]
        action_val = gm[2][end]
        @test all(isapprox.(path[2, :][(end - 5):end], 0, atol = 1e-3))
        @test all(isapprox.(action_val, 0.3375, atol = 1e-3))
    end

    @testset "HeymannVandenEijnden" begin # broken
        # method = "HeymannVandenEijnden"
        # gm = geometric_min_action_method(sys, x_i, x_f, maxiter = 10, method=method)

        # path = gm[1][end]
        # action_val = gm[2][end]
        # @test all(isapprox.(path[2, :][(end - 5):end], 0, atol = 1e-3))
        # @test all(isapprox.(action_val, 0.3375, atol = 1e-3))
    end # HeymannVandenEijnden
end # gMAM Meier Stein
"""
