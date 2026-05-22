using CriticalTransitions, StaticArrays
using Test
using LinearAlgebra
using Random

const CT = CriticalTransitions

@testset "DiagonalNoise update_p!: 1D OU multiplicative" begin
    α = 0.3
    b1d(u, p, t) = SA[-u[1]]
    g1d(u, p, t) = SA[sqrt(1 + α * u[1]^2);;]
    ds = CoupledSDEs(b1d, SA[1.0]; g = g1d, noise_prototype = SMatrix{1, 1}(0.0))

    sys = FreidlinWentzellHamiltonian(ds)
    @test sys isa FreidlinWentzellHamiltonian{<:Any, 1, <:Any, <:Any, <:Any, DiagonalNoise}

    Nt = 40
    x0 = reshape(collect(range(1.0, -1.0; length = Nt)), 1, Nt)
    s_arc = range(0; stop = 1, length = Nt)
    α_arc = zeros(Nt)
    CT.interpolate_path!(x0, α_arc, s_arc)
    xdot = zeros(size(x0))
    p = zeros(size(x0))
    λ = zeros(1, Nt)
    CT.central_diff!(xdot, x0)
    CT.update_p!(p, λ, x0, xdot, sys)
    @test all(isfinite, λ)
    @test all(isfinite, p)
end

@testset "GeneralNoise update_p!: 2D off-diagonal" begin
    b2(u, p, t) = SA[-u[1], -u[2]]
    function g2(u, p, t)
        s11 = 1 + 0.2 * u[1]; s22 = 1 + 0.2 * u[2]
        s12 = 0.3 * u[2];     s21 = 0.3 * u[1]
        return @SMatrix [s11 s12; s21 s22]
    end
    ds = CoupledSDEs(
        b2, SA[1.0, 0.0]; g = g2, noise_prototype = SMatrix{2, 2}(zeros(2, 2)),
    )
    sys = FreidlinWentzellHamiltonian(ds)
    @test sys isa FreidlinWentzellHamiltonian{<:Any, 2, <:Any, <:Any, <:Any, GeneralNoise}

    Nt = 40
    xx = collect(range(1.0, 0.0; length = Nt))
    yy = collect(range(0.0, 1.0; length = Nt))
    x0 = Matrix([xx yy]')
    s_arc = range(0; stop = 1, length = Nt)
    α_arc = zeros(Nt)
    CT.interpolate_path!(x0, α_arc, s_arc)
    xdot = zeros(size(x0))
    p = zeros(size(x0))
    λ = zeros(1, Nt)
    CT.central_diff!(xdot, x0)
    CT.update_p!(p, λ, x0, xdot, sys)
    @test all(isfinite, λ)
    @test all(isfinite, p)
end

@testset "sgMAM end-to-end: 1D OU multiplicative converges" begin
    α = 0.3
    b1d(u, p, t) = SA[-u[1]]
    g1d(u, p, t) = SA[sqrt(1 + α * u[1]^2);;]
    ds = CoupledSDEs(b1d, SA[1.0]; g = g1d, noise_prototype = SMatrix{1, 1}(0.0))
    sys = FreidlinWentzellHamiltonian(ds)

    Nt = 80
    x_initial = reshape(collect(range(1.0, -1.0; length = Nt)), 1, Nt)
    res = minimize_simple_geometric_action(
        sys, x_initial, GeometricGradient(; stepsize = 1.0);
        maxiters = 500, show_progress = false,
    )
    @test isfinite(res.action)
end

@testset "sgMAM end-to-end: 2D off-diagonal multiplicative converges" begin
    Random.seed!(0)
    b2(u, p, t) = SA[-u[1], -u[2]]
    function g2(u, p, t)
        s11 = 1 + 0.2 * u[1]; s22 = 1 + 0.2 * u[2]
        s12 = 0.3 * u[2];     s21 = 0.3 * u[1]
        return @SMatrix [s11 s12; s21 s22]
    end
    ds = CoupledSDEs(
        b2, SA[1.0, 0.0]; g = g2, noise_prototype = SMatrix{2, 2}(zeros(2, 2)),
    )
    sys = FreidlinWentzellHamiltonian(ds)

    Nt = 60
    xx = collect(range(1.0, 0.0; length = Nt))
    yy = collect(range(0.0, 1.0; length = Nt))
    x_initial = Matrix([xx yy]')
    res = minimize_simple_geometric_action(
        sys, x_initial, GeometricGradient(; stepsize = 1.0);
        maxiters = 300, show_progress = false,
    )
    @test isfinite(res.action)
end
