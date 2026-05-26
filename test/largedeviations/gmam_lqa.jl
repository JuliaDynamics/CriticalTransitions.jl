using Test
using CriticalTransitions
using LinearAlgebra
using StaticArrays

# Van der Pol from paper Sec 6.2.
function vdp_drift(u, p, t)
    x, y = u
    return SA[x - x^3 / 3 + y - y^3 / 9, x + 0.9]
end

function vdp_orbit(; Nτ = 400)
    sys0 = CoupledODEs(vdp_drift, [2.0, 0.0])
    pre = trajectory(sys0, 200.0; Δt = 0.01)[1]
    u0 = collect(pre[end])
    sys1 = CoupledODEs(vdp_drift, u0)
    T_period = 5.7966
    tr = trajectory(sys1, T_period - T_period / Nτ; Δt = T_period / Nτ)[1]
    pts = StateSpaceSet([SA[tr[k][1], tr[k][2]] for k in 1:Nτ])
    return pts, T_period
end

@testset "gMAM-LQA initial guess lies on tube surface" begin
    # The first node must be `γ(τ_start) + h · Ẽ(τ_start) · ẑ` with `‖ẑ‖ = 1`, so the
    # Cartesian distance to γ(τ_start) equals `h` exactly. `ẑ` should be the smallest
    # eigenvalue eigenvector of `G(τ_start)`.
    σ = 1.0
    pts, T = vdp_orbit(; Nτ = 200)
    sys = CoupledSDEs(vdp_drift, collect(pts[1]); noise_strength = σ)
    lc = LimitCycleFrame(pts, T, sys)
    G = local_quasipotential(lc)
    x_f = SA[2.0, -2.5]
    h = 0.05
    npoints = 50
    path0 = CriticalTransitions._gmam_lqa_initial_path(lc, G, x_f; tube_radius = h, npoints = npoints)
    @test size(path0) == (2, npoints)
    @test isapprox(path0[:, end], collect(x_f); atol = 1.0e-10)
    # Identify τ_start: orbit point closest to x_f.
    k0 = argmin([norm(collect(lc.γ[k]) - collect(x_f)) for k in 1:length(lc.γ)])
    γ_k0 = collect(lc.γ[k0])
    # Tube-surface invariant: ‖x₁ - γ(τ_start)‖ = h.
    @test isapprox(norm(path0[:, 1] - γ_k0), h; atol = 1.0e-12)
    # ẑ pulled back through Ẽ should match the smallest-eigenvalue eigenvector of G(k0).
    ẑ_back = lc.Ẽ[:, :, k0]' * (path0[:, 1] - γ_k0) / h
    @test isapprox(norm(ẑ_back), 1.0; atol = 1.0e-10)
    Gk0 = Symmetric((G[:, :, k0] + G[:, :, k0]') / 2)
    F = eigen(Gk0)
    ẑ_min = F.vectors[:, 1]
    @test isapprox(abs(dot(ẑ_back, ẑ_min)), 1.0; atol = 1.0e-10)
end

@testset "gMAM-LQA inner sweep reduces gMAM action" begin
    σ = 1.0
    pts, T = vdp_orbit(; Nτ = 200)
    sys = CoupledSDEs(vdp_drift, collect(pts[1]); noise_strength = σ)
    lc = LimitCycleFrame(pts, T, sys)
    G = local_quasipotential(lc)
    x_f = SA[2.0, -2.5]
    path = CriticalTransitions._gmam_lqa_initial_path(lc, G, x_f; tube_radius = 0.05, npoints = 80)
    ws = CriticalTransitions.geometric_gradient_workspace(sys, path)
    S0 = geometric_action(sys, path)
    sweep_ws = CriticalTransitions._InnerSweepWS(path)
    CriticalTransitions._gmam_lqa_inner_sweep!(
        path, ws, sys, sweep_ws, GeometricGradient(; stepsize = 0.01); maxiters = 50,
    )
    S1 = geometric_action(sys, path)
    # Observed drop is ~27% on this fixture; assert at least 20% as a guardrail against regressions.
    @test S1 < 0.8 * S0
end

@testset "Boundary step is monotone on the augmented action" begin
    # With a starting state offset from the optimal (k0, ẑ), one boundary-step iteration
    # must not increase the augmented action. (For d=2 the normal plane is 1D so we test
    # decrease from an offset k0 rather than a rotated ẑ.)
    σ = 1.0
    pts, T = vdp_orbit(; Nτ = 200)
    sys = CoupledSDEs(vdp_drift, collect(pts[1]); noise_strength = σ)
    lc = LimitCycleFrame(pts, T, sys)
    G = local_quasipotential(lc)
    x_f = SA[2.0, -2.5]
    k0_best = argmin([norm(collect(lc.γ[k]) - collect(x_f)) for k in 1:200])
    k0 = mod1(k0_best + 1, 200)  # one grid point off the closest match
    Gk = Symmetric((G[:, :, k0] + G[:, :, k0]') / 2)
    ẑ = eigen(Gk).vectors[:, 1]
    h = 0.05
    x1 = collect(lc.γ[k0]) + h * (lc.Ẽ[:, :, k0] * ẑ)
    npoints = 80
    path = zeros(2, npoints)
    for k in 1:npoints
        s = (k - 1) / (npoints - 1)
        path[:, k] = (1 - s) * x1 + s * collect(x_f)
    end
    state = (k0, ẑ)
    S0 = 0.5 * h^2 * dot(ẑ, Gk * ẑ) + geometric_action(sys, path)
    state2, path2 = CriticalTransitions._gmam_lqa_boundary_step(state, lc, G, sys, path, x_f, h)
    Gk2 = Symmetric((G[:, :, state2[1]] + G[:, :, state2[1]]') / 2)
    S1 = 0.5 * h^2 * dot(state2[2], Gk2 * state2[2]) + geometric_action(sys, path2)
    @test S1 ≤ S0
    # Step must move toward k0_best (or stay put): the new index is within {k0-1, k0, k0+1}.
    @test state2[1] ∈ (mod1(k0 - 1, 200), k0, mod1(k0 + 1, 200))
end

@testset "gMAM-LQA matches Stuart-Landau closed-form V(r)=(r²-μ)²/2" begin
    # Stuart-Landau has analytic quasi-potential V(r) = (r² - μ)²/2 in trace-normalized
    # action units. This validates the full gMAM-LQA pipeline against a closed form
    # (NOT just the leading-order LQA term ½·G·R², which is only asymptotic in R → 0).
    function stuart_landau_drift_test(u, p, t)
        x, y = u; μ_, ω_ = p; r2 = x^2 + y^2
        return SA[(μ_ - r2) * x - ω_ * y, (μ_ - r2) * y + ω_ * x]
    end
    μ, ω = 1.0, 1.0
    sys = CoupledSDEs(stuart_landau_drift_test, zeros(2), (μ, ω); noise_strength = 1.0)
    T_period = 2π / ω
    Nτ = 400
    τ_grid = range(0.0, T_period; length = Nτ + 1)[1:Nτ]
    pts = StateSpaceSet([SA[sqrt(μ) * cos(ω * t), sqrt(μ) * sin(ω * t)] for t in τ_grid])
    lc = LimitCycleFrame(pts, T_period, sys)
    G = local_quasipotential(lc)
    R = 0.3
    x_f = SA[sqrt(μ) + R, 0.0]
    V_exact = ((sqrt(μ) + R)^2 - μ)^2 / 2
    mp = minimize_geometric_action(
        sys, lc, x_f, GeometricGradient(; stepsize = 0.002);
        G = G, tube_radius = 0.02, npoints = 300,
        inner_maxiters = 500, maxiters = 400, reltol = 1.0e-9,
    )
    # gMAM-LQA computes the full nonlinear quasi-potential, not just the leading-order LQA.
    # Observed precision at these settings is ~0.6%; assert 1.5% as the regression guardrail.
    @test isapprox(mp.action, V_exact; rtol = 0.015)
end

@testset "gMAM-LQA reproduces paper Sec. 6.2 action ≈ 0.1567" begin
    # Cameron 2012 (Sec. 6.2 reference): minimum-action path from Γ to saddle (-0.9, 0.6942)
    # has quasi-potential 0.1567. Paper Lin et al. (2018) reports ≈0.1599 with N=160.
    σ = 1.0
    pts, T = vdp_orbit(; Nτ = 400)
    sys = CoupledSDEs(vdp_drift, collect(pts[1]); noise_strength = σ)
    lc = LimitCycleFrame(pts, T, sys)
    G = local_quasipotential(lc)
    x_saddle = SA[-0.9, 0.6942]
    mp = minimize_geometric_action(
        sys, lc, x_saddle, GeometricGradient(; stepsize = 0.005);
        G = G, tube_radius = 0.05, npoints = 160,
        inner_maxiters = 300, maxiters = 200, reltol = 1.0e-7,
    )
    # Observed: 0.1558 (relerr 0.6%). Tightened guardrail at 1.5%.
    @test isapprox(mp.action, 0.1567; rtol = 0.015)
end

@testset "gMAM-LQA with AdaptiveGeometricGradient matches GeometricGradient" begin
    σ = 1.0
    pts, T = vdp_orbit(; Nτ = 400)
    sys = CoupledSDEs(vdp_drift, collect(pts[1]); noise_strength = σ)
    lc = LimitCycleFrame(pts, T, sys)
    G = local_quasipotential(lc)
    x_saddle = SA[-0.9, 0.6942]
    mp_geom = minimize_geometric_action(
        sys, lc, x_saddle, GeometricGradient(; stepsize = 0.005);
        G = G, tube_radius = 0.05, npoints = 160,
        inner_maxiters = 300, maxiters = 200, reltol = 1.0e-7,
    )
    mp_adapt = minimize_geometric_action(
        sys, lc, x_saddle,
        AdaptiveGeometricGradient(; stepsize = 0.05, probe_length = 50);
        G = G, tube_radius = 0.05, npoints = 160,
        inner_maxiters = 300, maxiters = 200, reltol = 1.0e-7,
    )
    @test isapprox(mp_adapt.action, mp_geom.action; rtol = 0.02)
end
