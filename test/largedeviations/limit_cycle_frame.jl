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
    # Stuart-Landau closed-form multipliers: trivial 1, nontrivial exp(-2μT).
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
    # At γ(0)=(1,0), b = (-y, x) = (0, 1), so e₀ should equal (0, 1).
    @test isapprox(lc.E[:, 1, 1], [0.0, 1.0]; atol = 1.0e-6)
    # Normal direction at τ=0 is the radial direction (1, 0) (up to sign).
    @test isapprox(abs(lc.E[1, 2, 1]), 1.0; atol = 1.0e-5)
    @test isapprox(lc.E[2, 2, 1], 0.0; atol = 1.0e-5)
    # Periodicity (sign-permissive: orientable bundle, |e_j(T) - e_j(0)| within one Δτ step).
    Δθ = 2π / 200
    @test min(
        norm(lc.E[:, 2, 1] - lc.E[:, 2, end]),
        norm(lc.E[:, 2, 1] + lc.E[:, 2, end])
    ) < 2 * Δθ
end

@testset "Stuart-Landau M̃ ≈ -2μ and Ã ≈ 1 analytically" begin
    # SL radial linearization: ḋr = -2μ·dr ⇒ M̃ = -2μ.
    # Trace-normalised isotropic 2D noise gives a_normalized = I ⇒ Ã = 1 in the radial direction.
    μ, ω, σ = 1.0, 1.0, 0.7
    sys = CoupledSDEs(stuart_landau_drift, zeros(2), (μ, ω); noise_strength = σ)
    pts, T = stuart_landau_orbit(μ, ω; Nτ = 200)
    lc = LimitCycleFrame(pts, T, sys)
    @test all(isapprox.(lc.M̃[1, 1, :], -2 * μ; atol = 1.0e-4))
    @test all(isapprox.(lc.Ã[1, 1, :], 1.0; atol = 1.0e-9))
end

@testset "LimitCycleFrame anti-periodic detection" begin
    μ, ω = 1.0, 1.0
    sys = CoupledSDEs(stuart_landau_drift, zeros(2), (μ, ω); noise_strength = 0.1)
    pts, T = stuart_landau_orbit(μ, ω; Nτ = 50)
    @test_throws ArgumentError LimitCycleFrame(pts, T, sys; tol_periodic = 1.0e-16)
end
