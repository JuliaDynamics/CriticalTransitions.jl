using Test
using CriticalTransitions
using LinearAlgebra
using StaticArrays
using OrdinaryDiffEq: Tsit5

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
    k0 = argmin([norm(collect(lc.γ[k]) - collect(x_f)) for k in 1:length(lc.γ)])
    path0 = CriticalTransitions._gmam_lqa_straight_init(lc, G, x_f, k0, h, npoints)
    @test size(path0) == (2, npoints)
    @test isapprox(path0[:, end], collect(x_f); atol = 1.0e-10)
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
    k0 = argmin([norm(collect(lc.γ[k]) - collect(x_f)) for k in 1:length(lc.γ)])
    path = CriticalTransitions._gmam_lqa_straight_init(lc, G, x_f, k0, 0.05, 80)
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
    # For `d = 2` the normal plane is 1D so we test decrease from an offset `τ`
    # rather than a rotated `ẑ`.
    σ = 1.0
    pts, T = vdp_orbit(; Nτ = 200)
    sys = CoupledSDEs(vdp_drift, collect(pts[1]); noise_strength = σ)
    lc = LimitCycleFrame(pts, T, sys)
    G = local_quasipotential(lc)
    x_f = SA[2.0, -2.5]
    k0_best = argmin([norm(collect(lc.γ[k]) - collect(x_f)) for k in 1:200])
    k0 = mod1(k0_best + 1, 200)
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
    τ_0 = (k0 - 1) * T / 200
    state = (τ_0, ẑ)
    S0 = 0.5 * h^2 * dot(ẑ, Gk * ẑ) + geometric_action(sys, path)
    state2, path2 = CriticalTransitions._gmam_lqa_boundary_step_continuous(
        state, lc, G, sys, path, x_f, h,
    )
    Gτ2 = CriticalTransitions._G_at(G, state2[1], T, lc.S)
    S1 = 0.5 * h^2 * dot(state2[2], Gτ2 * state2[2]) + geometric_action(sys, path2)
    @test S1 ≤ S0
    @test mod(state2[1] - τ_0 + T / 2, T) - T / 2 ≤ T / 200 + 1.0e-12
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

@testset "gMAM-LQA initial-k0 scan beats the closest-γ heuristic on bimodal S(k0)" begin
    # 4D forced/dissipative slow-flow (coupled parametric RLC resonators in the
    # rotating frame): S(k0) is bimodal so the closest-γ heuristic lands in the
    # wrong basin (~1e-7) and the boundary walk cannot tunnel; the global scan
    # finds the true minimum (~3e-8).
    rgamma = 49268 / (3223137 * 2π)
    rj = -1435.0e9 / (3223137 * 2π)^2
    romega = 3216000 * 2π / (3223137 * 2π)
    ralpha = -1.0
    rlam = 2 * 0.0967 / ((3223137 * 2π) / 49268 * 0.052)
    κ = rj / (2 * romega)
    δ = romega / 2 - 1 / (2 * romega)
    μp = rlam * 1 / (4 * romega)
    cf = 3 * ralpha / (8 * romega)
    g = rgamma / 2
    function slowflow(u, p, t)
        u1, v1, u2, v2 = u[1], u[2], u[3], u[4]
        r1 = u1^2 + v1^2; r2 = u2^2 + v2^2
        return SA[
            κ * v2 - g * u1 + (δ - μp) * v1 - cf * r1 * v1,
            -κ * u2 - g * v1 - (δ + μp) * u1 + cf * r1 * u1,
            κ * v1 - g * u2 + δ * v2 - cf * r2 * v2,
            -κ * u1 - g * v2 - δ * u2 + cf * r2 * u2,
        ]
    end

    # Hard-coded point on Γ_+ and refined period from a one-off OptimizedShooting run
    # of the same drift (kept out of the test path to avoid the PeriodicOrbits / LM
    # dependency). The orbit-sampling step below traces Γ deterministically.
    u_on_lc = [
        -0.020069690676890243, 0.021455699931653267,
        0.01078726666608763, -0.03591127646718661,
    ]
    T_period = 4718.437300530521
    Nτ = 200
    # Tight Tsit5 because the period is ~4700; default tolerances accumulate enough
    # phase drift over a single revolution that the LC sample does not close.
    sys_det = CoupledODEs(
        slowflow, u_on_lc;
        diffeq = (alg = Tsit5(), abstol = 1.0e-12, reltol = 1.0e-12),
    )
    tr = trajectory(sys_det, T_period - T_period / Nτ; Δt = T_period / Nτ)[1]
    pts = StateSpaceSet([SA[tr[k][1], tr[k][2], tr[k][3], tr[k][4]] for k in 1:Nτ])
    sys = CoupledSDEs(slowflow, u_on_lc; noise_strength = 1.0)
    lc = LimitCycleFrame(pts, T_period, sys)
    G = local_quasipotential(lc; maxiters = 1500)

    x_f = zeros(4)
    opt = AdaptiveGeometricGradient(;
        stepsize = 1.0e4, stepsize_min = 1.0e-12, stepsize_max = 1.0e6,
        shrink = 0.5, grow = 1.1, probe_length = 400,
    )

    # `scan_maxiters = 500` is the minimum we have validated; below ~300 the
    # bimodality has not yet surfaced in the post-sweep action.
    mp_auto = minimize_geometric_action(
        sys, lc, x_f, opt;
        G = G, tube_radius = 0.002, npoints = 200,
        inner_maxiters = 500, maxiters = 100, reltol = 1.0e-9,
        scan_maxiters = 500,
    )
    @test mp_auto.action < 6.0e-8

    γs = [collect(pts[k]) for k in 1:Nτ]
    k0_legacy = argmin(k -> norm(γs[k] - x_f), 1:Nτ)
    mp_legacy = minimize_geometric_action(
        sys, lc, x_f, opt;
        G = G, tube_radius = 0.002, npoints = 200,
        inner_maxiters = 500, maxiters = 100, reltol = 1.0e-9,
        k0_init = k0_legacy,
    )
    @test mp_legacy.action > 9.0e-8

    @test mp_auto.action < mp_legacy.action
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
