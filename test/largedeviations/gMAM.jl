using CriticalTransitions
using Test

function _fhn_gmam(u, p, t)
    x, y = u
    ϵ, β, α, γ, κ, I = p
    return SA[(-α * x^3 + γ * x - κ * y + I) / ϵ, -β * y + x]
end

@testset "gMAM FitzHugh-Nagumo" begin
    p = [0.1, 3, 1, 1, 1, 0]
    σ = 0.1
    fhn = CoupledSDEs(_fhn_gmam, zeros(2), p; noise_strength = σ)
    x_i = SA[sqrt(2 / 3), sqrt(2 / 27)]
    x_f = SA[0.001, 0.0]
    res = minimize_geometric_action(
        fhn, x_i, x_f; npoints = 30, maxiters = 500, show_progress = false
    )
    S = geometric_action(fhn, Matrix(res.path)')
    @test isapprox(S, 0.18, atol = 0.01)
end

@testset "GeometricGradient" begin
    function meier_stein(u, p, t)
        x, y = u
        dx = x - x^3 - 10 * x * y^2
        dy = -(1 + x^2) * y
        return SA[dx, dy]
    end
    σ = 0.25
    sys = CoupledSDEs(meier_stein, zeros(2); noise_strength = σ)

    xx = range(-1.0, 1.0; length = 30)
    yy = 0.3 .* (-xx .^ 2 .+ 1)
    init = Matrix([xx yy]')

    gm = minimize_geometric_action(
        sys, init, GeometricGradient(); maxiters = 500, verbose = false, show_progress = false
    )

    path = Matrix(gm.path)'
    action_val = gm.action
    @test all(isapprox.(path[2, :][(end - 5):end], 0, atol = 1.0e-3))
    @test all(isapprox.(action_val, 0.3375, atol = 1.0e-3))
end

@testset "gMAM backtracking" begin
    function meier_stein(u, p, t)
        x, y = u
        dx = x - x^3 - 10 * x * y^2
        dy = -(1 + x^2) * y
        return SA[dx, dy]
    end
    σ = 0.25
    sys = CoupledSDEs(meier_stein, zeros(2); noise_strength = σ)

    xx = range(-1.0, 1.0; length = 30)
    yy = 0.3 .* (-xx .^ 2 .+ 1)
    init = Matrix([xx yy]')

    res_bt = minimize_geometric_action(
        sys, init, GeometricGradient(; stepsize = 1.0e6); maxiters = 100, show_progress = false
    )
    @test isfinite(res_bt.action)

    actions = Float64[]
    for ss in [0.01, 1.0, 1.0e3]
        res = minimize_geometric_action(
            sys, init, GeometricGradient(; stepsize = ss); maxiters = 500, show_progress = false
        )
        push!(actions, res.action)
    end
    @test maximum(actions) / minimum(actions) < 1.05
end

@testset "GeometricGradient constructor" begin
    opt = GeometricGradient()
    @test opt.stepsize isa Float64
    @test opt.max_backtracks == 10

    opt2 = GeometricGradient(; stepsize = 1, shrink = 0.5)
    @test opt2.stepsize isa Float64
    @test opt2.stepsize == 1.0

    opt3 = GeometricGradient(; max_backtracks = 0, stepsize = 42.0)
    @test opt3.max_backtracks == 0
    @test opt3.stepsize == 42.0
end
