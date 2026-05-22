using CriticalTransitions
using ModelingToolkitBase
using Test
using LinearAlgebra
using StaticArrays

const CT = CriticalTransitions

@testset "FreidlinWentzellHamiltonian inner-loop type-stable" begin
    function meier_stein_ts(u, p, t)
        x, y = u
        return SA[x - x^3 - 10 * x * y^2, -(1 + x^2) * y]
    end
    ds = CoupledSDEs(meier_stein_ts, zeros(2); noise_strength = 0.25)
    sys = FreidlinWentzellHamiltonian(ds)

    Nt = 20
    xx = range(-1.0, 1.0; length = Nt)
    yy = 0.3 .* (-xx .^ 2 .+ 1)
    xpath = Matrix([xx yy]')
    xdot = zeros(size(xpath)); ppath = zeros(size(xpath)); λ = zeros(1, Nt)
    pdot = zeros(size(xpath)); xdotdot = zeros(size(xpath))
    CT.central_diff!(xdot, xpath)
    CT.update_p!(ppath, λ, xpath, xdot, sys)
    Hx = sys.H_x(xpath, ppath)
    CT.central_diff!(pdot, ppath); CT.central_diff!(xdotdot, xdot)

    @inferred CT.update_p!(ppath, λ, xpath, xdot, sys)
    @inferred CT.update_x!(xpath, λ, pdot, xdotdot, Hx, sys, 1.0)
end

@testset "AdditiveNoise update_x! no per-iteration sparse alloc" begin
    function meier_stein(u, p, t)
        x, y = u
        return SA[x - x^3 - 10 * x * y^2, -(1 + x^2) * y]
    end
    ds = CoupledSDEs(meier_stein, zeros(2); noise_strength = 0.25)
    sys = FreidlinWentzellHamiltonian(ds)
    Nt = 60
    xx = range(-1.0, 1.0; length = Nt)
    yy = 0.3 .* (-xx .^ 2 .+ 1)
    x_initial = Matrix([xx yy]')
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

@testset "FreidlinWentzellHamiltonian carries NoiseShape" begin
    f_lin(u, p, t) = SA[-u[1], -u[2]]
    ds_ode = CoupledODEs(f_lin, SA[0.0, 0.0])
    sys_ode = FreidlinWentzellHamiltonian(ds_ode)
    @test sys_ode isa FreidlinWentzellHamiltonian{<:Any, 2, <:Any, <:Any, <:Any, AdditiveNoise}
    @test sys_ode.a(zeros(2)) ≈ LinearAlgebra.Diagonal(ones(2))

    ds_iso = CoupledSDEs(f_lin, SA[0.0, 0.0]; noise_strength = 1.0)
    sys_iso = FreidlinWentzellHamiltonian(ds_iso)
    @test sys_iso isa FreidlinWentzellHamiltonian{<:Any, 2, <:Any, <:Any, <:Any, AdditiveNoise}

    H_x_user(x, p) = zeros(size(x))
    H_p_user(x, p) = ones(size(x))
    sys_user = FreidlinWentzellHamiltonian{false, 2}(H_x_user, H_p_user)
    @test sys_user isa FreidlinWentzellHamiltonian{false, 2, <:Any, <:Any, <:Any, AdditiveNoise}
end

@testset "FreidlinWentzellHamiltonian KPO" begin
    λ = 3 / 1.21 * 2 / 295
    ω0 = 1.0
    ω = 1.0
    γ = 1 / 295
    η = 0
    α = -1

    function fu(u, v)
        return (-4 * γ * ω * u - 2 * λ * v - 4 * (ω0 - ω^2) * v - 3 * α * v * (u^2 + v^2)) /
            (8 * ω)
    end
    function fv(u, v)
        return (-4 * γ * ω * v - 2 * λ * u + 4 * (ω0 - ω^2) * u + 3 * α * u * (u^2 + v^2)) /
            (8 * ω)
    end

    dfvdv(u, v) = (-4 * γ * ω + 6 * α * u * v) / (8 * ω)
    dfudu(u, v) = (-4 * γ * ω - 6 * α * u * v) / (8 * ω)
    dfvdu(u, v) = (-2 * λ + 4 * (ω0 - ω^2) + 9 * α * u^2 + 3 * α * v^2) / (8 * ω)
    dfudv(u, v) = (-2 * λ - 4 * (ω0 - ω^2) - 3 * α * u^2 - 9 * α * v^2) / (8 * ω)

    function H_x(x, p) # ℜ² → ℜ²
        u, v = eachrow(x)
        pu, pv = eachrow(p)

        H_u = @. pu * dfudu(u, v) + pv * dfvdu(u, v)
        H_v = @. pu * dfudv(u, v) + pv * dfvdv(u, v)
        return Matrix([H_u H_v]')
    end
    function H_p(x, p) # ℜ² → ℜ²
        u, v = eachrow(x)
        pu, pv = eachrow(p)

        H_pu = @. pu + fu(u, v)
        H_pv = @. pv + fv(u, v)
        return Matrix([H_pu H_pv]')
    end

    @independent_variables t
    D = Differential(t)
    sts = @variables u(t) v(t)

    eqs = [D(u) ~ fu(u, v), D(v) ~ fv(u, v)]
    @mtkcompile sysMTK = System(eqs, t)
    prob = ODEProblem(sysMTK, Dict(sts .=> zeros(2)), (0.0, 100.0); jac = true)
    ds = CoupledODEs(prob)

    sys = FreidlinWentzellHamiltonian{false, 2}(H_x, H_p)
    sys′ = FreidlinWentzellHamiltonian(ds)

    Nt = 500  # number of discrete time steps
    p_r = rand(2, Nt)
    x_r = rand(2, Nt)

    @test sys′.H_x(x_r, p_r) ≈ sys.H_x(x_r, p_r)
    @test sys′.H_p(x_r, p_r) ≈ sys.H_p(x_r, p_r)
end

@testset "FreidlinWentzellHamiltonian MTK" begin
    @independent_variables t
    D = Differential(t)
    sts = @variables u(t) v(t)

    @parameters λ = 3 / 1.21 * 2 / 295 ω0 = 1.0 ω = 1.0 γ = 1 / 295 η = 0 α = -1

    eqs = [
        D(u) ~
            (-4 * γ * ω * u - 2 * λ * v - 4 * (ω0 - ω^2) * v - 3 * α * v * (u^2 + v^2)) /
            (8 * ω),
        D(v) ~
            (-4 * γ * ω * v - 2 * λ * u + 4 * (ω0 - ω^2) * u + 3 * α * u * (u^2 + v^2)) /
            (8 * ω),
    ]
    @mtkcompile sysMTK = System(eqs, t)
    prob = ODEProblem(sysMTK, Dict(sts .=> zeros(2)), (0.0, 100.0); jac = true)
    ds = CoupledODEs(prob)
    sys = FreidlinWentzellHamiltonian(ds)

    @test sys.H_x(zeros(2), zeros(2)) ≈ zeros(2)
    @test sys.H_p(zeros(2), zeros(2)) ≈ zeros(2)
end

@testset "sgMAM GeometricGradient" begin
    H_x(x, p) = zeros(size(x))
    H_p(x, p) = ones(size(x))
    sys = FreidlinWentzellHamiltonian{false, 2}(H_x, H_p)

    xx = collect(range(-1.0, 1.0; length = 20))
    yy = 0.3 .* (-xx .^ 2 .+ 1)
    x_initial = Matrix([xx yy]')

    res_small = minimize_geometric_action(
        sys,
        x_initial,
        GeometricGradient(; stepsize = 1.0e-6, max_backtracks = 0);
        maxiters = 2,
        show_progress = false,
    )
    res_large = minimize_geometric_action(
        sys,
        x_initial,
        GeometricGradient(; stepsize = 1.0, max_backtracks = 0);
        maxiters = 2,
        show_progress = false,
    )

    @test res_small.action != res_large.action
end

@testset "sgMAM GeometricGradient Backtracking" begin
    CT = CriticalTransitions

    function meier_stein(u, p, t) # out-of-place
        x, y = u
        dx = x - x^3 - 10 * x * y^2
        dy = -(1 + x^2) * y
        return SA[dx, dy]
    end

    σ = 0.25
    ds = CoupledSDEs(meier_stein, zeros(2); noise_strength = σ)
    sys = FreidlinWentzellHamiltonian(ds)

    xx = range(-1.0, 1.0; length = 60)
    yy = 0.3 .* (-xx .^ 2 .+ 1)
    x_initial = Matrix([xx yy]')

    # Baseline action on an arclength-reparameterized path
    x0 = deepcopy(x_initial)
    Nt = size(x0, 2)
    s = range(0; stop = 1, length = Nt)
    α = zeros(Nt)
    CT.interpolate_path!(x0, α, s)
    xdot = zeros(size(x0))
    p = zeros(size(x0))
    λ = zeros(1, Nt)
    CT.central_diff!(xdot, x0)
    CT.update_p!(p, λ, x0, xdot, sys)
    S0 = CT.FW_action(xdot, p)
    @test isfinite(S0)

    # Single iteration with huge stepsize: backtracking should prevent blowup
    opt = GeometricGradient(; max_backtracks = 20, stepsize = 1.0e6)
    res = minimize_geometric_action(
        sys, x_initial, opt; maxiters = 1, show_progress = false
    )
    @test isfinite(res.action)
    @test res.action <= S0 + 1.0e-10

    # Multi-iteration: backtracking should converge to a lower action than the initial path
    bt_max_backtracks = 20
    bt_maxiters = 200
    res_bt = minimize_geometric_action(
        sys,
        x_initial,
        GeometricGradient(; max_backtracks = bt_max_backtracks, stepsize = 1.0);
        maxiters = bt_maxiters,
        show_progress = false,
    )
    @test res_bt.action < S0

    # Fairness baseline: no-backtracking gets the same total trial-step budget
    # as the backtracking run could consume (max_backtracks+1 trials per outer iter).
    res_no_bt = minimize_geometric_action(
        sys,
        x_initial,
        GeometricGradient(; max_backtracks = 0, stepsize = 1.0);
        maxiters = bt_maxiters * (bt_max_backtracks + 1),
        show_progress = false,
    )
    # Backtracking should reach at least as good (or better) action
    @test res_bt.action <= res_no_bt.action + 1.0e-6
end

@testset "sgMAM action matches geometric_action on converged path" begin
    # Regression guard for the spurious /2 in `FW_action`: once converged, the
    # path-action reported by sgMAM and the geometric action evaluated on the
    # same path must agree (both are S_FW on the instanton). Before the fix,
    # sgMAM reported half the geometric action.
    function meier_stein(u, p, t)
        x, y = u
        dx = x - x^3 - 10 * x * y^2
        dy = -(1 + x^2) * y
        return SA[dx, dy]
    end
    σ = 0.25
    ds = CoupledSDEs(meier_stein, zeros(2); noise_strength = σ)
    sys = FreidlinWentzellHamiltonian(ds)

    xx = range(-1.0, 1.0; length = 60)
    yy = 0.3 .* (-xx .^ 2 .+ 1)
    x_initial = Matrix([xx yy]')

    res = minimize_geometric_action(
        sys,
        x_initial,
        GeometricGradient(; max_backtracks = 20, stepsize = 1.0);
        maxiters = 500,
        show_progress = false,
    )

    S_geo = geometric_action(ds, Matrix(res.path)', 1.0)
    @test isapprox(res.action, S_geo; rtol = 1.0e-3)
end

@testset "sgMAM step-size insensitivity" begin
    function meier_stein(u, p, t)
        x, y = u
        dx = x - x^3 - 10 * x * y^2
        dy = -(1 + x^2) * y
        return SA[dx, dy]
    end
    σ = 0.25
    ds = CoupledSDEs(meier_stein, zeros(2); noise_strength = σ)
    sys = FreidlinWentzellHamiltonian(ds)

    xx = range(-1.0, 1.0; length = 60)
    yy = 0.3 .* (-xx .^ 2 .+ 1)
    x_initial = Matrix([xx yy]')

    # Different starting stepsizes should converge to the same action
    actions = Float64[]
    for ss in [1.0, 100.0, 1.0e4]
        res = minimize_geometric_action(
            sys,
            x_initial,
            GeometricGradient(; stepsize = ss);
            maxiters = 500,
            show_progress = false,
        )
        push!(actions, res.action)
    end
    # All actions should be within 5% of each other
    @test maximum(actions) / minimum(actions) < 1.05
end

@testset "sgMAM convergence tolerances" begin
    function meier_stein(u, p, t)
        x, y = u
        dx = x - x^3 - 10 * x * y^2
        dy = -(1 + x^2) * y
        return SA[dx, dy]
    end
    σ = 0.25
    ds = CoupledSDEs(meier_stein, zeros(2); noise_strength = σ)
    sys = FreidlinWentzellHamiltonian(ds)

    xx = range(-1.0, 1.0; length = 60)
    yy = 0.3 .* (-xx .^ 2 .+ 1)
    x_initial = Matrix([xx yy]')

    # With a tight reltol, should converge and get a finite action
    res_tol = minimize_geometric_action(
        sys,
        x_initial,
        GeometricGradient(; stepsize = 100.0);
        maxiters = 10_000,
        reltol = 1.0e-6,
        show_progress = false,
        verbose = false,
    )
    @test isfinite(res_tol.action)

    # Without reltol and more iterations should get same or better action
    res_notol = minimize_geometric_action(
        sys,
        x_initial,
        GeometricGradient(; stepsize = 100.0);
        maxiters = 10_000,
        show_progress = false,
    )
    @test res_notol.action <= res_tol.action + 1.0e-10
end

function _maier_stein_setup(; Nt = 60)
    function meier_stein(u, p, t)
        x, y = u
        dx = x - x^3 - 10 * x * y^2
        dy = -(1 + x^2) * y
        return SA[dx, dy]
    end
    ds = CoupledSDEs(meier_stein, zeros(2); noise_strength = 0.25)
    sys = FreidlinWentzellHamiltonian(ds)
    xx = range(-1.0, 1.0; length = Nt)
    yy = 0.3 .* (-xx .^ 2 .+ 1)
    return ds, sys, Matrix([xx yy]')
end

@testset "AdaptiveGeometricGradient constructor" begin
    opt = AdaptiveGeometricGradient()
    @test opt isa AdaptiveGeometricGradient
    @test opt isa CriticalTransitions.GMAMOptimizer
    @test opt.probe_length == 200
    @test 0 < opt.shrink < 1
    @test opt.grow > 1

    opt2 = AdaptiveGeometricGradient(;
        stepsize = 50.0, probe_length = 100, shrink = 0.4, grow = 1.5,
    )
    @test opt2.stepsize == 50.0
    @test opt2.probe_length == 100
    @test opt2.shrink == 0.4
    @test opt2.grow == 1.5

    # Validation
    @test_throws ArgumentError AdaptiveGeometricGradient(; probe_length = 0)
    @test_throws ArgumentError AdaptiveGeometricGradient(; shrink = 0.0)
    @test_throws ArgumentError AdaptiveGeometricGradient(; shrink = 1.0)
    @test_throws ArgumentError AdaptiveGeometricGradient(; grow = 1.0)
end

@testset "AdaptiveGeometricGradient Maier–Stein" begin
    CT = CriticalTransitions
    _, sys, x_initial = _maier_stein_setup()

    # Baseline action
    x0 = deepcopy(x_initial)
    Nt = size(x0, 2)
    s = range(0; stop = 1, length = Nt)
    α = zeros(Nt)
    CT.interpolate_path!(x0, α, s)
    xdot = zeros(size(x0)); p = zeros(size(x0)); λ = zeros(1, Nt)
    CT.central_diff!(xdot, x0)
    CT.update_p!(p, λ, x0, xdot, sys)
    S0 = CT.FW_action(xdot, p)
    @test isfinite(S0)

    # Adaptive should strictly improve on the initial path
    res_ad = minimize_geometric_action(
        sys, x_initial,
        AdaptiveGeometricGradient(; stepsize = 100.0, probe_length = 50);
        maxiters = 500, show_progress = false,
    )
    @test isfinite(res_ad.action)
    @test res_ad.action < S0

    # Adaptive should be competitive with backtracking GeometricGradient
    res_gg = minimize_geometric_action(
        sys, x_initial,
        GeometricGradient(; stepsize = 100.0);
        maxiters = 1000, show_progress = false,
    )
    # Within 10% (Maier–Stein is overdamped: both regimes find similar paths;
    # adaptive uses 2× compute per probe window so with the same maxiters it may
    # be slightly behind).
    @test res_ad.action <= res_gg.action * 1.1 + 1.0e-8
end

@testset "AdaptiveGeometricGradient stepsize insensitivity" begin
    _, sys, x_initial = _maier_stein_setup()

    actions = Float64[]
    for ss in (1.0, 100.0, 1.0e4)
        res = minimize_geometric_action(
            sys, x_initial,
            AdaptiveGeometricGradient(; stepsize = ss, probe_length = 50);
            maxiters = 500, show_progress = false,
        )
        @test isfinite(res.action)
        push!(actions, res.action)
    end
    # All actions within 5% — the adaptive controller should erase the initial-step bias
    @test maximum(actions) / minimum(actions) < 1.05
end

@testset "AdaptiveGeometricGradient convergence tolerances" begin
    _, sys, x_initial = _maier_stein_setup()

    # reltol-based early stop
    res_tol = minimize_geometric_action(
        sys, x_initial,
        AdaptiveGeometricGradient(; stepsize = 100.0, probe_length = 50);
        maxiters = 10_000, reltol = 1.0e-6, show_progress = false,
    )
    @test isfinite(res_tol.action)

    res_notol = minimize_geometric_action(
        sys, x_initial,
        AdaptiveGeometricGradient(; stepsize = 100.0, probe_length = 50);
        maxiters = 10_000, show_progress = false,
    )
    @test res_notol.action <= res_tol.action + 1.0e-8
end

@testset "AdaptiveGeometricGradient API forms" begin
    _, sys, x_initial = _maier_stein_setup(; Nt = 40)

    # Accepts a StateSpaceSet
    sss = StateSpaceSet(x_initial')
    res_sss = minimize_geometric_action(
        sys, sss,
        AdaptiveGeometricGradient(; stepsize = 100.0, probe_length = 30);
        maxiters = 100, show_progress = false,
    )
    @test isfinite(res_sss.action)
end

@testset "AdaptiveGeometricGradient underdamped KPO" begin
    # KPO in the underdamped regime is the canonical test for AdaptiveGeometricGradient:
    # backtracking GeometricGradient settles at stepsize_max which over-smooths the
    # spiral structure near the attractor, while the adaptive variant detects this over
    # the probe window and shrinks the step.
    λ_val = 3 / 1.21 * 2 / 100
    α_val = -1.0
    γ_val = λ_val / 2 * 0.05  # κ = 0.05, strongly underdamped
    ω0 = 1.0; ω_v = 1.0

    fu(u, v) = (-4γ_val * ω_v * u - 2λ_val * v - 4(ω0 - ω_v^2) * v - 3α_val * v * (u^2 + v^2)) / (8ω_v)
    fv(u, v) = (-4γ_val * ω_v * v - 2λ_val * u + 4(ω0 - ω_v^2) * u + 3α_val * u * (u^2 + v^2)) / (8ω_v)
    dfudu(u, v) = (-4γ_val * ω_v - 6α_val * u * v) / (8ω_v)
    dfudv(u, v) = (-2λ_val - 4(ω0 - ω_v^2) - 3α_val * u^2 - 9α_val * v^2) / (8ω_v)
    dfvdu(u, v) = (-2λ_val + 4(ω0 - ω_v^2) + 9α_val * u^2 + 3α_val * v^2) / (8ω_v)
    dfvdv(u, v) = (-4γ_val * ω_v + 6α_val * u * v) / (8ω_v)

    function H_x(x, p)
        u, v = eachrow(x); pu, pv = eachrow(p)
        H_u = @. pu * dfudu(u, v) + pv * dfvdu(u, v)
        H_v = @. pu * dfudv(u, v) + pv * dfvdv(u, v)
        return Matrix([H_u H_v]')
    end
    function H_p(x, p)
        u, v = eachrow(x); pu, pv = eachrow(p)
        H_pu = @. pu + fu(u, v)
        H_pv = @. pv + fv(u, v)
        return Matrix([H_pu H_pv]')
    end
    sys = FreidlinWentzellHamiltonian{false, 2}(H_x, H_p)

    κ = 2γ_val / λ_val
    r = sqrt(2λ_val * sqrt(1 - κ^2) / (3 * abs(α_val)))
    θ = atan(-κ, -sqrt(1 - κ^2)) / 2
    xa = [r * cos(θ), r * sin(θ)]; xb = -xa
    Nt = 200
    s = collect(range(0; stop = 1, length = Nt))
    xx = @. (xb[1] - xa[1]) * s + xa[1]
    yy = @. (xb[2] - xa[2]) * s + xa[2] + 0.01sin(2π * s)
    x_initial = Matrix([xx yy]')

    # Backtracking GeometricGradient (incumbent) starting from a large step size
    res_gg = minimize_geometric_action(
        sys, x_initial,
        GeometricGradient(; stepsize = 100.0);
        maxiters = 2000, show_progress = false,
    )

    # Adaptive — should reach a lower or equal action on this underdamped problem
    res_ad = minimize_geometric_action(
        sys, x_initial,
        AdaptiveGeometricGradient(; stepsize = 100.0, probe_length = 100);
        maxiters = 2000, show_progress = false,
    )
    @test isfinite(res_ad.action)
    @test isfinite(res_gg.action)
    # On underdamped systems, adaptive should be strictly competitive (no worse than 0.5%)
    @test res_ad.action <= res_gg.action * 1.005 + 1.0e-8

    # Returned path is well-formed
    @test length(res_ad.path) == Nt
    @test size(res_ad.generalized_momentum) == (2, Nt)
    @test size(res_ad.path_velocity) == (2, Nt)
    @test size(res_ad.λ) == (1, Nt)
end
