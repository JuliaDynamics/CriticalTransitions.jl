using CriticalTransitions, StaticArrays
using Test
using LinearAlgebra
using Random

const CT = CriticalTransitions

# Canonical 1D OU with multiplicative diagonal noise σ(x) = √(1 + α x²).
_b1d(u, p, t) = SA[-u[1]]
_make_1d_ou(α) = let
    g(u, p, t) = SA[sqrt(1 + α * u[1]^2);;]
    CoupledSDEs(_b1d, SA[1.0]; g = g, noise_prototype = SMatrix{1, 1}(0.0))
end

# Canonical 2D linear drift with full off-diagonal multiplicative noise.
_b2(u, p, t) = SA[-u[1], -u[2]]
function _g2(u, p, t)
    s11 = 1 + 0.2 * u[1]; s22 = 1 + 0.2 * u[2]
    s12 = 0.3 * u[2];     s21 = 0.3 * u[1]
    return @SMatrix [s11 s12; s21 s22]
end
_make_2d_offdiag() =
    CoupledSDEs(_b2, SA[1.0, 0.0]; g = _g2, noise_prototype = SMatrix{2, 2}(zeros(2, 2)))

@testset "DiagonalNoise update_p!: 1D OU multiplicative" begin
    ds = _make_1d_ou(0.3)

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
    ds = _make_2d_offdiag()
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
    sys = FreidlinWentzellHamiltonian(_make_1d_ou(0.3))
    Nt = 80
    x_initial = reshape(collect(range(1.0, -1.0; length = Nt)), 1, Nt)
    res = minimize_geometric_action(
        sys, x_initial, GeometricGradient(; stepsize = 1.0);
        maxiters = 500, show_progress = false,
    )
    @test isfinite(res.action)
end

@testset "sgMAM end-to-end: 2D off-diagonal multiplicative converges" begin
    Random.seed!(0)
    sys = FreidlinWentzellHamiltonian(_make_2d_offdiag())
    Nt = 60
    xx = collect(range(1.0, 0.0; length = Nt))
    yy = collect(range(0.0, 1.0; length = Nt))
    x_initial = Matrix([xx yy]')
    res = minimize_geometric_action(
        sys, x_initial, GeometricGradient(; stepsize = 1.0);
        maxiters = 300, show_progress = false,
    )
    @test isfinite(res.action)
end

@testset "_action_metric: additive returns constant inv" begin
    f_lin(u, p, t) = SA[-u[1], -u[2]]
    ds = CoupledSDEs(f_lin, SA[0.0, 0.0]; noise_strength = 2.0)
    metric = CT._action_metric(ds)
    @test metric isa Base.Returns
    @test metric(zeros(2)) ≈ inv(LinearAlgebra.Diagonal([1.0, 1.0]))
end

@testset "_action_metric: state-dependent evaluates per point" begin
    metric = CT._action_metric(_make_1d_ou(0.3))
    @test !(metric isa Base.Returns)
    @test metric([1.0])[1, 1] ≈ 1.0
end

@testset "fw_action: 1D OU multiplicative vs analytic Simpson" begin
    α = 0.3
    ds = _make_1d_ou(α)

    N, T = 200, 1.0
    path = reduce(hcat, range([1.0], [0.0]; length = N))
    time = range(0.0, T; length = N)
    S = fw_action(ds, path, time)

    function simpson(f, a, b, n)
        h = (b - a) / n; s = f(a) + f(b)
        for i in 1:2:(n - 1); s += 4 * f(a + i * h); end
        for i in 2:2:(n - 2); s += 2 * f(a + i * h); end
        return s * h / 3
    end
    s_norm = 1 + α
    integrand = t -> t^2 / (1 + α * (1 - t)^2)
    analytic = s_norm * simpson(integrand, 0, 1, 10_000) / 2
    @test isapprox(S, analytic; rtol = 1.0e-4)
end

@testset "gMAM diagonal multiplicative converges" begin
    Random.seed!(0)
    ds = _make_1d_ou(0.3)
    Nt = 80
    x_initial = reshape(collect(range(1.0, -1.0; length = Nt)), 1, Nt)
    res_g = minimize_geometric_action(
        ds, x_initial, GeometricGradient(; stepsize = 1.0);
        maxiters = 500, show_progress = false,
    )
    @test isfinite(res_g.action)
end

@testset "Additive non-diagonal Σ goes through GeneralNoise" begin
    function meier_stein(u, p, t)
        x, y = u
        return SA[x - x^3 - 10 * x * y^2, -(1 + x^2) * y]
    end
    Nt = 60
    xx = range(-1.0, 1.0; length = Nt)
    yy = 0.3 .* (-xx .^ 2 .+ 1)
    x_initial = Matrix([xx yy]')

    θ = 0.4
    R = [cos(θ) -sin(θ); sin(θ) cos(θ)]
    D = Diagonal([0.5, 2.0])
    Q_rot = R * D * R'
    ds_rot = CoupledSDEs(meier_stein, zeros(2); covariance = Q_rot)
    sys_rot = FreidlinWentzellHamiltonian(ds_rot)
    @test sys_rot isa FreidlinWentzellHamiltonian{<:Any, 2, <:Any, <:Any, <:Any, GeneralNoise}

    res = minimize_geometric_action(
        sys_rot, x_initial, GeometricGradient(; stepsize = 1.0);
        maxiters = 500, show_progress = false,
    )
    @test isfinite(res.action)
end

@testset "gMAM general multiplicative converges" begin
    Random.seed!(0)
    ds = _make_2d_offdiag()
    Nt = 60
    xx = collect(range(1.0, 0.0; length = Nt))
    yy = collect(range(0.0, 1.0; length = Nt))
    x_initial = Matrix([xx yy]')
    res_g = minimize_geometric_action(
        ds, x_initial, GeometricGradient(; stepsize = 1.0);
        maxiters = 300, show_progress = false,
    )
    @test isfinite(res_g.action)
end
