using CriticalTransitions, StaticArrays
using Test
using LinearAlgebra
using Random

function _make_1d_ou(α)
    b(u, p, t) = SA[-u[1]]
    g(u, p, t) = SA[sqrt(1 + α * u[1]^2);;]
    return CoupledSDEs(
        b, SA[1.0]; g = g, noise_prototype = SMatrix{1, 1}(0.0),
    )
end

function _make_2d_offdiag()
    b(u, p, t) = SA[-u[1], -u[2]]
    function g(u, p, t)
        s11 = 1 + 0.2 * u[1]
        s22 = 1 + 0.2 * u[2]
        s12 = 0.3 * u[2]
        s21 = 0.3 * u[1]
        return SMatrix{2, 2}(s11, s21, s12, s22)
    end
    return CoupledSDEs(
        b, SA[1.0, 0.0]; g = g,
        noise_prototype = SMatrix{2, 2}(zeros(2, 2)),
    )
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

@testset "fw_action: 1D OU multiplicative vs analytic Simpson" begin
    α = 0.3
    ds = _make_1d_ou(α)

    N, T = 200, 1.0
    path = reduce(hcat, range([1.0], [0.0]; length = N))
    time = range(0.0, T; length = N)
    S = fw_action(ds, path, time)

    function simpson(f, a, b, n)
        h = (b - a) / n
        s = f(a) + f(b)
        for i in 1:2:(n - 1)
            s += 4 * f(a + i * h)
        end
        for i in 2:2:(n - 2)
            s += 2 * f(a + i * h)
        end
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

@testset "MAM (FW) on multiplicative noise agrees with gMAM at large T" begin
    Random.seed!(0)
    ds = _make_1d_ou(0.3)
    Nt = 80

    xinit_g = reshape(collect(range(1.0, -1.0; length = Nt)), 1, Nt)
    res_g = minimize_geometric_action(
        ds, xinit_g, GeometricGradient(; stepsize = 1.0);
        maxiters = 2000, show_progress = false,
    )
    S_g = res_g.action

    init = reduce(hcat, range([1.0], [-1.0]; length = Nt))
    Ts = [0.5, 1.0, 4.0, 16.0]
    S_mam = map(Ts) do T
        res = minimize_action(ds, init, T; maxiters = 2000, show_progress = false)
        return res.action
    end

    @test issorted(S_mam; rev = true)
    @test S_mam[1] > 3 * S_g
    @test S_mam[end] ≈ S_g rtol = 0.05
    @test_throws ArgumentError minimize_action(
        ds, init, 1.0;
        functional = "OM", noise_strength = 0.1,
        maxiters = 5, show_progress = false,
    )
end

@testset "sgMAM additive non-diagonal covariance converges" begin
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

@testset "gMAM state-dep a ≡ I matches additive a ≡ I" begin
    Random.seed!(0)
    function ms_drift(u, p, t)
        x, y = u
        return SA[x - x^3 - 10 * x * y^2, -(1 + x^2) * y]
    end

    ds_add = CoupledSDEs(ms_drift, zeros(2); noise_strength = 1.0)

    function g_rotation(u, p, t)
        c = cos(0.5 * u[1])
        s = sin(0.5 * u[1])
        return @SMatrix [c -s; s c]
    end
    ds_sd = CoupledSDEs(
        ms_drift, zeros(2); g = g_rotation,
        noise_prototype = SMatrix{2, 2}(zeros(2, 2)),
    )

    Nt = 40
    xx = range(-1.0, 1.0; length = Nt)
    yy = 0.3 .* (-xx .^ 2 .+ 1)
    x_initial = Matrix([xx yy]')

    res_add = minimize_geometric_action(
        ds_add, x_initial, GeometricGradient(; stepsize = 1.0);
        maxiters = 500, show_progress = false,
    )
    res_sd = minimize_geometric_action(
        ds_sd, x_initial, GeometricGradient(; stepsize = 1.0);
        maxiters = 500, show_progress = false,
    )

    @test isapprox(res_add.action, res_sd.action; rtol = 1.0e-6)
end
