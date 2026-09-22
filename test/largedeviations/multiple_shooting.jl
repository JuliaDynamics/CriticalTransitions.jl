using CriticalTransitions
using Test
using LinearAlgebra

using Logging
using Statistics: mean

bistable_rule(u, p, t) = SA[u[1] - u[1]^3]
maier_stein_rule(u, p, t) = SA[u[1] - u[1]^3 - u[1] * u[2]^2, -(1 + u[1]^2) * u[2]]

bistable_H() = FreidlinWentzellHamiltonian(CoupledSDEs(bistable_rule, [0.0]; noise_strength = 1.0))
maier_stein_H() = FreidlinWentzellHamiltonian(CoupledSDEs(maier_stein_rule, [0.0, 0.0]; noise_strength = 1.0))

linear_ramp(x_start, x_end, Nt) = Matrix(reshape(collect(range(x_start, x_end; length = Nt)), 1, Nt))

maier_stein_arc(x_start, x_end, Nt; amp = 0.2) = let
    xs = collect(range(x_start, x_end; length = Nt))
    ys = [amp * sin(π * (x - x_start) / (x_end - x_start)) for x in xs]
    Matrix(vcat(xs', ys'))
end

run_silent(args...; kwargs...) = with_logger(NullLogger()) do
    minimize_geometric_action(args...; kwargs..., show_progress = false)
end

@testset "MultipleShooting constructor" begin
    opt = MultipleShooting()
    @test opt isa MultipleShooting
    @test opt.nshoots == 10
    opt2 = MultipleShooting(; nshoots = 5, abstol = 1.0e-7, reltol = 1.0e-5)
    @test opt2.nshoots == 5
    @test opt2.abstol == 1.0e-7
    @test_throws ArgumentError MultipleShooting(; nshoots = 1)
end

@testset "1D bistable additive: attractor → saddle (half instanton)" begin
    H = bistable_H()
    Nt = 20
    x_init = linear_ramp(-1.0, 0.0, Nt)

    res = run_silent(H, x_init, MultipleShooting(; nshoots = 6, maxiters = 200))

    cfg = [res.path[k][1] for k in 1:length(res.path)]
    ramp = collect(range(-1.0, 0.0; length = Nt))
    @test maximum(abs.(cfg .- ramp)) < 5.0e-2
    @test res.path[1][1] ≈ -1.0 atol = 1.0e-6
    @test abs(res.path[end][1]) < 5.0e-2
    @test isapprox(res.λ, 1.0; rtol = 1.0e-2)
end

@testset "1D bistable: full instanton via splitting (no interior fp)" begin
    H = bistable_H()
    Nt = 20

    res_l = run_silent(H, linear_ramp(-1.0, 0.0, Nt), MultipleShooting(; nshoots = 6, maxiters = 200))
    # +1 → 0 is the noise-driven half (mirror of −1 → 0); saddle → attractor is deterministic.
    res_r = run_silent(H, linear_ramp(1.0, 0.0, Nt), MultipleShooting(; nshoots = 6, maxiters = 200))

    @test isapprox(res_l.action, res_r.action; rtol = 1.0e-6)
end

@testset "Convergence in nshoots" begin
    H = bistable_H()
    x_init = linear_ramp(-1.0, 0.0, 30)
    actions = map([4, 6, 10]) do n
        run_silent(H, x_init, MultipleShooting(; nshoots = n, maxiters = 200)).action
    end
    @test (maximum(actions) - minimum(actions)) / mean(actions) < 5.0e-3
end

@testset "Maier-Stein additive: shooting matches sgMAM (attractor → saddle)" begin
    H_ms = maier_stein_H()
    x_init = maier_stein_arc(-1.0, 0.0, 30)
    res_shoot = run_silent(H_ms, x_init, MultipleShooting(; nshoots = 8, maxiters = 200, abstol = 1.0e-6))
    res_sg = run_silent(H_ms, x_init, GeometricGradient(); maxiters = 5000)
    @test isapprox(res_shoot.action, res_sg.action; rtol = 1.0e-2)
end

@testset "1D bistable multiplicative (DiagonalNoise)" begin
    γ = 0.5
    bistable_diff(u, p, t) = SA[sqrt(1.0 + γ * u[1]^2);;]
    ds = CoupledSDEs(bistable_rule, [0.0]; g = bistable_diff, noise_prototype = zeros(1, 1))
    H = FreidlinWentzellHamiltonian(ds)
    x_init = linear_ramp(-1.0, 0.0, 30)

    run_silent(H, x_init, MultipleShooting(; nshoots = 8, maxiters = 200, abstol = 1.0e-6))
end

@testset "GeneralNoise: Maier-Stein with rotated constant Σ" begin
    θ = π / 6
    Q_rot = let c = cos(θ), s = sin(θ)
        R = [c -s; s c]
        R * Diagonal([1.0, 2.0]) * R'
    end
    ds = CoupledSDEs(maier_stein_rule, [0.0, 0.0]; covariance = Q_rot)
    H_rot = FreidlinWentzellHamiltonian(ds)
    x_init = maier_stein_arc(-1.0, 0.0, 30)
    res_shoot = run_silent(H_rot, x_init, MultipleShooting(; nshoots = 8, maxiters = 200, abstol = 1.0e-6))
    res_sg = run_silent(H_rot, x_init, GeometricGradient(); maxiters = 5000)
    @test isapprox(res_shoot.action, res_sg.action; rtol = 5.0e-2)
end

@testset "Generalized momentum on-shell quality" begin
    H = bistable_H()
    Nt = 30
    res = run_silent(H, linear_ramp(-1.0, 0.0, Nt), MultipleShooting(; nshoots = 6, maxiters = 200))
    @test size(res.generalized_momentum) == (Nt, 1)
end

@testset "Maier-Stein left-right symmetry" begin
    H_ms = maier_stein_H()
    res_l = run_silent(H_ms, maier_stein_arc(-1.0, 0.0, 30), MultipleShooting(; nshoots = 8, maxiters = 200, abstol = 1.0e-6))
    res_r = run_silent(H_ms, maier_stein_arc(1.0, 0.0, 30), MultipleShooting(; nshoots = 8, maxiters = 200, abstol = 1.0e-6))
    @test isapprox(res_l.action, res_r.action; rtol = 1.0e-3)
end
