using Test
using CriticalTransitions
using LinearAlgebra
using StaticArrays

# Stuart-Landau oscillator: stable limit cycle r = sqrt(μ), period 2π/ω.
# Analytic structure used as the reference everywhere in this file.
function stuart_landau_drift(u, p, t)
    x, y = u
    μ, ω = p
    r2 = x^2 + y^2
    return SA[(μ - r2) * x - ω * y, (μ - r2) * y + ω * x]
end

function stuart_landau_orbit(μ, ω; Nτ = 200)
    T = 2π / ω
    τ = range(0.0, T; length = Nτ + 1)[1:(end - 1)]
    r = sqrt(μ)
    pts = [SA[r * cos(ω * t), r * sin(ω * t)] for t in τ]
    return StateSpaceSet(pts), T
end

@testset "_state_transition_matrices initial value and identity at τ=0" begin
    μ, ω = 1.0, 1.0
    sys = CoupledSDEs(stuart_landau_drift, zeros(2), (μ, ω); noise_strength = 0.1)
    pts, T = stuart_landau_orbit(μ, ω; Nτ = 200)
    Φs = CriticalTransitions._state_transition_matrices(sys, pts, T)
    @test size(Φs) == (2, 2, 200)
    @test Φs[:, :, 1] ≈ I(2) atol = 1.0e-10
end

@testset "Stuart-Landau monodromy has Floquet multipliers (1, exp(-2μT))" begin
    μ, ω = 1.0, 1.0
    T = 2π / ω
    sys = CoupledSDEs(stuart_landau_drift, zeros(2), (μ, ω); noise_strength = 0.1)
    pts, _ = stuart_landau_orbit(μ, ω; Nτ = 400)
    tands = CriticalTransitions.DynamicalSystemsBase.TangentDynamicalSystem(
        CoupledODEs(sys); u0 = collect(pts[1]),
    )
    CriticalTransitions.DynamicalSystemsBase.step!(tands, T, true)
    Φ_M = Matrix(CriticalTransitions.DynamicalSystemsBase.current_deviations(tands))
    mults = sort(abs.(eigvals(Φ_M)))
    @test isapprox(mults[2], 1.0; atol = 1.0e-5)
    @test isapprox(mults[1], exp(-2 * μ * T); rtol = 1.0e-4)
end

@testset "LimitCycleFrame on Stuart-Landau: analytic frame at τ=0" begin
    μ, ω = 1.0, 1.0
    sys = CoupledSDEs(stuart_landau_drift, zeros(2), (μ, ω); noise_strength = 0.1)
    pts, T = stuart_landau_orbit(μ, ω; Nτ = 200)
    lc = LimitCycleFrame(pts, T, sys)
    # Purely periodic bundle: F = (+1,), sign mask S = (+1,).
    @test lc.F == (Int8(1),)
    @test lc.S == SMatrix{1, 1, Int8}(1)
    # At γ(0)=(1,0), b = (-y, x) = (0, 1), so e₀ should equal (0, 1).
    @test isapprox(lc.E[:, 1, 1], [0.0, 1.0]; atol = 1.0e-6)
    # Normal direction at τ=0 is the radial direction (1, 0) (up to sign).
    @test isapprox(abs(lc.E[1, 2, 1]), 1.0; atol = 1.0e-5)
    @test isapprox(lc.E[2, 2, 1], 0.0; atol = 1.0e-5)
    # Periodicity: with F=+I and smear, Ẽ closes T-periodically up to one Δτ step.
    Δθ = 2π / 200
    @test min(
        norm(lc.E[:, 2, 1] - lc.E[:, 2, end]),
        norm(lc.E[:, 2, 1] + lc.E[:, 2, end])
    ) < 2 * Δθ
end

@testset "Stuart-Landau M̃ ≈ -2μ and Ã ≈ 1 analytically" begin
    μ, ω, σ = 1.0, 1.0, 0.7
    sys = CoupledSDEs(stuart_landau_drift, zeros(2), (μ, ω); noise_strength = σ)
    pts, T = stuart_landau_orbit(μ, ω; Nτ = 200)
    lc = LimitCycleFrame(pts, T, sys)
    @test all(isapprox.(lc.M̃[1, 1, :], -2 * μ; atol = 1.0e-4))
    @test all(isapprox.(lc.Ã[1, 1, :], 1.0; atol = 1.0e-9))
end

@testset "LimitCycleFrame on 4D system with rotation-by-π in transverse block" begin
    # 4D flow: Stuart-Landau cycle in (x, y) plus a skew (z₁, z₂) block whose monodromy
    # equals exp(-αT)·R(π) = -exp(-αT)·I (rotation by π, orientable, det = +1).
    # Previously this was rejected by `@test_throws "anti-periodic"`. With the Schur-based
    # algorithm the rotation is absorbed into the smear; F = (+1, +1, +1) and the Riccati
    # produces the analytic answer G = diag(1, 1, 4): two z-directions with G = 2α = 1,
    # one radial direction with G = 4μ = 4.
    function skew_flow(u, p, t)
        x, y, z1, z2 = u
        μ, ω, α = p
        r2 = x^2 + y^2
        return SA[
            (μ - r2) * x - ω * y,
            (μ - r2) * y + ω * x,
            -α * z1 - (ω / 2) * z2,
            (ω / 2) * z1 - α * z2,
        ]
    end
    μ, ω, α = 1.0, 1.0, 0.5
    sys = CoupledSDEs(
        skew_flow, [sqrt(μ), 0.0, 0.0, 0.0], (μ, ω, α); noise_strength = 0.1,
    )
    T = 2π / ω
    Nτ = 200
    τ_grid = range(0.0, T; length = Nτ + 1)[1:Nτ]
    pts = StateSpaceSet(
        [
            SA[sqrt(μ) * cos(ω * t), sqrt(μ) * sin(ω * t), 0.0, 0.0] for t in τ_grid
        ]
    )
    lc = LimitCycleFrame(pts, T, sys)
    @test lc.F == (Int8(1), Int8(1), Int8(1))
    @test lc.S == ones(Int8, 3, 3)
    G = local_quasipotential(lc; maxiters = 300)
    @test all(isposdef(Symmetric(G[:, :, k])) for k in 1:Nτ)
    # Stationary G in the canonical basis: two transverse z-directions give G = 2α,
    # radial direction gives G = 4μ. Ordering depends on Schur, so sort eigenvalues.
    eigs = sort(eigvals(Symmetric(G[:, :, 1])))
    @test isapprox(eigs[1], 2 * α; atol = 1.0e-3)
    @test isapprox(eigs[2], 2 * α; atol = 1.0e-3)
    @test isapprox(eigs[3], 4 * μ; atol = 1.0e-3)
end

@testset "local_quasipotential two-leg path with forced antiperiodic F" begin
    # Genuine non-orientable transverse bundles cannot arise in autonomous ODEs on Rᵈ
    # (the bundle is always orientable). To exercise the two-leg integration path, we
    # construct a Stuart-Landau LC and forcibly set F = (-1,) on the struct. Since
    # F² = 1 in the m=1 case, the sign-mask S is still +1 and the two-leg path must
    # converge to the same fixed point as the standard single-leg integration.
    μ, ω, σ = 1.0, 1.0, 0.7
    sys = CoupledSDEs(stuart_landau_drift, zeros(2), (μ, ω); noise_strength = σ)
    pts, T = stuart_landau_orbit(μ, ω; Nτ = 200)
    lc = LimitCycleFrame(pts, T, sys)
    lc_anti = typeof(lc)(
        lc.γ, lc.period, lc.E, lc.Ẽ, lc.M̃, lc.Ã,
        (Int8(-1),), SMatrix{1, 1, Int8}(1),
    )
    @test any(lc_anti.F .== -1)
    G_anti = local_quasipotential(lc_anti)
    G_std = local_quasipotential(lc)
    @test maximum(abs.(G_anti .- G_std)) < 1.0e-6
end

@testset "G_at recovers the second half via the sign mask" begin
    μ, ω, σ = 1.0, 1.0, 0.7
    sys = CoupledSDEs(stuart_landau_drift, zeros(2), (μ, ω); noise_strength = σ)
    pts, T = stuart_landau_orbit(μ, ω; Nτ = 200)
    lc = LimitCycleFrame(pts, T, sys)
    G = local_quasipotential(lc)
    # F = (+1,) so the second half equals the first.
    @test G_at(lc, G, 0.5) ≈ G_at(lc, G, T + 0.5)
    @test G_at(lc, G, 0.5) ≈ G_at(lc, G, 2 * T + 0.5)
end

@testset "Warm-start G0 converges to the same fixed point" begin
    μ, ω, σ = 1.0, 1.0, 0.7
    sys = CoupledSDEs(stuart_landau_drift, zeros(2), (μ, ω); noise_strength = σ)
    pts, T = stuart_landau_orbit(μ, ω; Nτ = 200)
    lc = LimitCycleFrame(pts, T, sys)
    G_cold = local_quasipotential(lc)
    G_warm = local_quasipotential(lc; G0 = G_cold)
    @test maximum(abs.(G_warm .- G_cold)) < 1.0e-9
end
