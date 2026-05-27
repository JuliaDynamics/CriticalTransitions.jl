using Test
using CriticalTransitions
using LinearAlgebra
using StaticArrays

# Stuart-Landau fixture comes from limit_cycle_frame.jl included earlier.
@testset "PRDE RHS vanishes at analytic fixed-point on Stuart-Landau" begin
    μ, ω, σ = 1.0, 1.0, 0.7
    sys = CoupledSDEs(stuart_landau_drift, zeros(2), (μ, ω); noise_strength = σ)
    pts, T = stuart_landau_orbit(μ, ω; Nτ = 200)
    lc = LimitCycleFrame(pts, T, sys)
    Ã0 = lc.Ã[1, 1, 1]
    G_expected = 4 * μ / Ã0
    Gmat = fill(G_expected, 1, 1)
    out = similar(Gmat)
    CriticalTransitions._prde_rhs!(out, Gmat, lc, 0.0)
    @test isapprox(out[1, 1], 0.0; atol = 1.0e-4)
end

@testset "local_quasipotential converges to Stuart-Landau analytic G" begin
    μ, ω, σ = 1.0, 1.0, 0.7
    sys = CoupledSDEs(stuart_landau_drift, zeros(2), (μ, ω); noise_strength = σ)
    pts, T = stuart_landau_orbit(μ, ω; Nτ = 200)
    lc = LimitCycleFrame(pts, T, sys)
    G = local_quasipotential(lc)
    @test size(G) == (1, 1, 200)
    Ã0 = lc.Ã[1, 1, 1]
    G_expected = 4 * μ / Ã0
    @test all(isapprox.(G[1, 1, :], G_expected; rtol = 1.0e-5))
    for k in 1:200
        @test isposdef(Symmetric(G[:, :, k]))
    end
    @test isapprox(G[:, :, 1], G[:, :, end]; rtol = 1.0e-6)
end
