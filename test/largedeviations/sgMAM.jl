using CriticalTransitions
using Test
using LinearAlgebra

function _maier_stein_setup(; Nt = 60)
    function meier_stein(u, p, t)
        x, y = u
        return SA[x - x^3 - 10 * x * y^2, -(1 + x^2) * y]
    end
    ds = CoupledSDEs(meier_stein, zeros(2); noise_strength = 0.25)
    sys = FreidlinWentzellHamiltonian(ds)
    xx = range(-1.0, 1.0; length = Nt)
    yy = 0.3 .* (-xx .^ 2 .+ 1)
    return ds, sys, Matrix([xx yy]')
end

@testset "sgMAM public-call allocations" begin
    _, sys, x_initial = _maier_stein_setup()
    minimize_geometric_action(
        sys, x_initial, GeometricGradient(; stepsize = 1.0);
        maxiters = 1, show_progress = false,
    )
    bytes_before = Base.gc_num().total_allocd
    minimize_geometric_action(
        sys, x_initial, GeometricGradient(; stepsize = 1.0);
        maxiters = 10, show_progress = false,
    )
    bytes_after = Base.gc_num().total_allocd
    @test (bytes_after - bytes_before) < 5_000_000
end

@testset "FreidlinWentzellHamiltonian public fields" begin
    f_lin(u, p, t) = SA[-u[1], -u[2]]
    ds_ode = CoupledODEs(f_lin, SA[0.0, 0.0])
    sys_ode = FreidlinWentzellHamiltonian(ds_ode)
    @test sys_ode isa FreidlinWentzellHamiltonian{<:Any, 2}
    @test sys_ode.a isa Base.Returns
    @test sys_ode.a(zeros(2)) ≈ LinearAlgebra.Diagonal(ones(2))

    ds_iso = CoupledSDEs(f_lin, SA[0.0, 0.0]; noise_strength = 1.0)
    sys_iso = FreidlinWentzellHamiltonian(ds_iso)
    @test sys_iso isa FreidlinWentzellHamiltonian{<:Any, 2}
    @test sys_iso.a isa Base.Returns

    H_x_user(x, p) = zeros(size(x))
    H_p_user(x, p) = ones(size(x))
    sys_user = FreidlinWentzellHamiltonian{false, 2}(H_x_user, H_p_user)
    @test sys_user isa FreidlinWentzellHamiltonian{false, 2}
    @test sys_user.a isa Base.Returns
end

@testset "FreidlinWentzellHamiltonian automatic derivatives" begin
    λ = 3 / 1.21 * 2 / 295
    ω0 = 1.0
    ω = 1.0
    γ = 1 / 295
    α = -1

    fu(u, v) = (-4 * γ * ω * u - 2 * λ * v - 4 * (ω0 - ω^2) * v - 3 * α * v * (u^2 + v^2)) / (8 * ω)
    fv(u, v) = (-4 * γ * ω * v - 2 * λ * u + 4 * (ω0 - ω^2) * u + 3 * α * u * (u^2 + v^2)) / (8 * ω)
    dfvdv(u, v) = (-4 * γ * ω + 6 * α * u * v) / (8 * ω)
    dfudu(u, v) = (-4 * γ * ω - 6 * α * u * v) / (8 * ω)
    dfvdu(u, v) = (-2 * λ + 4 * (ω0 - ω^2) + 9 * α * u^2 + 3 * α * v^2) / (8 * ω)
    dfudv(u, v) = (-2 * λ - 4 * (ω0 - ω^2) - 3 * α * u^2 - 9 * α * v^2) / (8 * ω)

    function H_x(x, p)
        u, v = eachrow(x)
        pu, pv = eachrow(p)
        H_u = @. pu * dfudu(u, v) + pv * dfvdu(u, v)
        H_v = @. pu * dfudv(u, v) + pv * dfvdv(u, v)
        return Matrix([H_u H_v]')
    end
    function H_p(x, p)
        u, v = eachrow(x)
        pu, pv = eachrow(p)
        H_pu = @. pu + fu(u, v)
        H_pv = @. pv + fv(u, v)
        return Matrix([H_pu H_pv]')
    end

    kpo_rhs(u, p, t) = SA[fu(u[1], u[2]), fv(u[1], u[2])]
    ds = CoupledODEs(kpo_rhs, zeros(2))
    sys = FreidlinWentzellHamiltonian{false, 2}(H_x, H_p)
    sys_auto = FreidlinWentzellHamiltonian(ds)

    p_r = rand(2, 500)
    x_r = rand(2, 500)
    @test sys_auto.H_x(x_r, p_r) ≈ sys.H_x(x_r, p_r)
    @test sys_auto.H_p(x_r, p_r) ≈ sys.H_p(x_r, p_r)
end

@testset "GeometricGradient" begin
    H_x(x, p) = zeros(size(x))
    H_p(x, p) = ones(size(x))
    sys = FreidlinWentzellHamiltonian{false, 2}(H_x, H_p)

    xx = collect(range(-1.0, 1.0; length = 20))
    yy = 0.3 .* (-xx .^ 2 .+ 1)
    x_initial = Matrix([xx yy]')

    res_small = minimize_geometric_action(
        sys, x_initial, GeometricGradient(; stepsize = 1.0e-6, max_backtracks = 0);
        maxiters = 2, show_progress = false,
    )
    res_large = minimize_geometric_action(
        sys, x_initial, GeometricGradient(; stepsize = 1.0, max_backtracks = 0);
        maxiters = 2, show_progress = false,
    )
    @test res_small.action != res_large.action
end

@testset "GeometricGradient backtracking" begin
    _, sys, x_initial = _maier_stein_setup()
    S0 = minimize_geometric_action(
        sys, x_initial, GeometricGradient(; stepsize = 1.0);
        maxiters = 0, show_progress = false,
    ).action
    @test isfinite(S0)

    res = minimize_geometric_action(
        sys, x_initial, GeometricGradient(; max_backtracks = 20, stepsize = 1.0e6);
        maxiters = 1, show_progress = false,
    )
    @test isfinite(res.action)
    @test res.action <= S0 + 1.0e-10

    res_bt = minimize_geometric_action(
        sys, x_initial, GeometricGradient(; max_backtracks = 20, stepsize = 1.0);
        maxiters = 200, show_progress = false,
    )
    res_no_bt = minimize_geometric_action(
        sys, x_initial, GeometricGradient(; max_backtracks = 0, stepsize = 1.0);
        maxiters = 4200, show_progress = false,
    )
    @test res_bt.action < S0
    @test res_bt.action <= res_no_bt.action + 1.0e-6
end

@testset "sgMAM action matches geometric_action" begin
    ds, sys, x_initial = _maier_stein_setup()
    res = minimize_geometric_action(
        sys, x_initial, GeometricGradient(; max_backtracks = 20, stepsize = 1.0);
        maxiters = 500, show_progress = false,
    )
    S_geo = geometric_action(ds, Matrix(res.path)', 1.0)
    @test isapprox(res.action, S_geo; rtol = 1.0e-3)
end

@testset "GeometricGradient step-size insensitivity" begin
    _, sys, x_initial = _maier_stein_setup()
    actions = Float64[]
    for ss in (1.0, 100.0, 1.0e4)
        res = minimize_geometric_action(
            sys, x_initial, GeometricGradient(; stepsize = ss);
            maxiters = 500, show_progress = false,
        )
        push!(actions, res.action)
    end
    @test maximum(actions) / minimum(actions) < 1.05
end

@testset "GeometricGradient convergence tolerances" begin
    _, sys, x_initial = _maier_stein_setup()
    res_tol = minimize_geometric_action(
        sys, x_initial, GeometricGradient(; stepsize = 100.0);
        maxiters = 10_000, reltol = 1.0e-6, show_progress = false, verbose = false,
    )
    res_notol = minimize_geometric_action(
        sys, x_initial, GeometricGradient(; stepsize = 100.0);
        maxiters = 10_000, show_progress = false,
    )
    @test isfinite(res_tol.action)
    @test res_notol.action <= res_tol.action + 1.0e-10
end

@testset "AdaptiveGeometricGradient constructor" begin
    opt = AdaptiveGeometricGradient()
    @test opt isa AdaptiveGeometricGradient
    @test opt.probe_length == 200
    @test 0 < opt.shrink < 1
    @test opt.grow > 1

    opt2 = AdaptiveGeometricGradient(; stepsize = 50.0, probe_length = 100, shrink = 0.4, grow = 1.5)
    @test opt2.stepsize == 50.0
    @test opt2.probe_length == 100
    @test opt2.shrink == 0.4
    @test opt2.grow == 1.5

    @test_throws ArgumentError AdaptiveGeometricGradient(; probe_length = 0)
    @test_throws ArgumentError AdaptiveGeometricGradient(; shrink = 0.0)
    @test_throws ArgumentError AdaptiveGeometricGradient(; shrink = 1.0)
    @test_throws ArgumentError AdaptiveGeometricGradient(; grow = 1.0)
end

@testset "AdaptiveGeometricGradient Maier-Stein" begin
    _, sys, x_initial = _maier_stein_setup()
    S0 = minimize_geometric_action(
        sys, x_initial, GeometricGradient(; stepsize = 1.0);
        maxiters = 0, show_progress = false,
    ).action

    res_ad = minimize_geometric_action(
        sys, x_initial, AdaptiveGeometricGradient(; stepsize = 100.0, probe_length = 50);
        maxiters = 500, show_progress = false,
    )
    res_gg = minimize_geometric_action(
        sys, x_initial, GeometricGradient(; stepsize = 100.0);
        maxiters = 1000, show_progress = false,
    )
    @test isfinite(res_ad.action)
    @test res_ad.action < S0
    @test res_ad.action <= res_gg.action * 1.1 + 1.0e-8
end

@testset "AdaptiveGeometricGradient step-size insensitivity" begin
    _, sys, x_initial = _maier_stein_setup()
    actions = Float64[]
    for ss in (1.0, 100.0, 1.0e4)
        res = minimize_geometric_action(
            sys, x_initial, AdaptiveGeometricGradient(; stepsize = ss, probe_length = 50);
            maxiters = 500, show_progress = false,
        )
        @test isfinite(res.action)
        push!(actions, res.action)
    end
    @test maximum(actions) / minimum(actions) < 1.05
end

@testset "AdaptiveGeometricGradient convergence tolerances" begin
    _, sys, x_initial = _maier_stein_setup()
    res_tol = minimize_geometric_action(
        sys, x_initial, AdaptiveGeometricGradient(; stepsize = 100.0, probe_length = 50);
        maxiters = 10_000, reltol = 1.0e-6, show_progress = false,
    )
    res_notol = minimize_geometric_action(
        sys, x_initial, AdaptiveGeometricGradient(; stepsize = 100.0, probe_length = 50);
        maxiters = 10_000, show_progress = false,
    )
    @test isfinite(res_tol.action)
    @test res_notol.action <= res_tol.action + 1.0e-8
end

@testset "AdaptiveGeometricGradient StateSpaceSet" begin
    _, sys, x_initial = _maier_stein_setup(; Nt = 40)
    res = minimize_geometric_action(
        sys, StateSpaceSet(x_initial'), AdaptiveGeometricGradient(; stepsize = 100.0, probe_length = 30);
        maxiters = 100, show_progress = false,
    )
    @test isfinite(res.action)
end

@testset "AdaptiveGeometricGradient underdamped KPO" begin
    λ_val = 3 / 1.21 * 2 / 100
    α_val = -1.0
    γ_val = λ_val / 2 * 0.05
    ω0 = 1.0
    ω_v = 1.0

    fu(u, v) = (-4γ_val * ω_v * u - 2λ_val * v - 4(ω0 - ω_v^2) * v - 3α_val * v * (u^2 + v^2)) / (8ω_v)
    fv(u, v) = (-4γ_val * ω_v * v - 2λ_val * u + 4(ω0 - ω_v^2) * u + 3α_val * u * (u^2 + v^2)) / (8ω_v)
    dfudu(u, v) = (-4γ_val * ω_v - 6α_val * u * v) / (8ω_v)
    dfudv(u, v) = (-2λ_val - 4(ω0 - ω_v^2) - 3α_val * u^2 - 9α_val * v^2) / (8ω_v)
    dfvdu(u, v) = (-2λ_val + 4(ω0 - ω_v^2) + 9α_val * u^2 + 3α_val * v^2) / (8ω_v)
    dfvdv(u, v) = (-4γ_val * ω_v + 6α_val * u * v) / (8ω_v)

    function H_x(x, p)
        u, v = eachrow(x)
        pu, pv = eachrow(p)
        H_u = @. pu * dfudu(u, v) + pv * dfvdu(u, v)
        H_v = @. pu * dfudv(u, v) + pv * dfvdv(u, v)
        return Matrix([H_u H_v]')
    end
    function H_p(x, p)
        u, v = eachrow(x)
        pu, pv = eachrow(p)
        H_pu = @. pu + fu(u, v)
        H_pv = @. pv + fv(u, v)
        return Matrix([H_pu H_pv]')
    end
    sys = FreidlinWentzellHamiltonian{false, 2}(H_x, H_p)

    κ = 2γ_val / λ_val
    r = sqrt(2λ_val * sqrt(1 - κ^2) / (3 * abs(α_val)))
    θ = atan(-κ, -sqrt(1 - κ^2)) / 2
    xa = [r * cos(θ), r * sin(θ)]
    xb = -xa
    Nt = 200
    s = collect(range(0; stop = 1, length = Nt))
    xx = @. (xb[1] - xa[1]) * s + xa[1]
    yy = @. (xb[2] - xa[2]) * s + xa[2] + 0.01sin(2π * s)
    x_initial = Matrix([xx yy]')

    res_gg = minimize_geometric_action(
        sys, x_initial, GeometricGradient(; stepsize = 100.0);
        maxiters = 2000, show_progress = false,
    )
    res_ad = minimize_geometric_action(
        sys, x_initial, AdaptiveGeometricGradient(; stepsize = 100.0, probe_length = 100);
        maxiters = 2000, show_progress = false,
    )
    @test isfinite(res_ad.action)
    @test isfinite(res_gg.action)
    @test res_ad.action <= res_gg.action * 1.005 + 1.0e-8
    @test length(res_ad.path) == Nt
    @test size(res_ad.generalized_momentum) == (2, Nt)
    @test size(res_ad.path_velocity) == (2, Nt)
    @test size(res_ad.λ) == (1, Nt)
end
