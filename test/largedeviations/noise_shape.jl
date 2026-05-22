using CriticalTransitions
using Test
using LinearAlgebra
using StaticArrays

@testset "NoiseShape hierarchy" begin
    @test NoiseShape === CriticalTransitions.NoiseShape
    @test AdditiveNoise <: NoiseShape
    @test DiagonalNoise <: NoiseShape
    @test GeneralNoise <: NoiseShape
    @test AdditiveNoise() isa NoiseShape
    @test DiagonalNoise() isa NoiseShape
    @test GeneralNoise() isa NoiseShape
end

@testset "_is_rank_deficient" begin
    @test !CriticalTransitions._is_rank_deficient(Matrix{Float64}(LinearAlgebra.I(3)))
    @test !CriticalTransitions._is_rank_deficient(LinearAlgebra.Diagonal([1.0, 2.0, 3.0]))
    @test  CriticalTransitions._is_rank_deficient(LinearAlgebra.Diagonal([1.0, 0.0, 3.0]))
    @test  CriticalTransitions._is_rank_deficient([1.0 1.0; 1.0 1.0])
end

@testset "_probe_points" begin
    u₀ = [0.5, -0.3]
    probes = CriticalTransitions._probe_points(u₀)
    @test length(probes) == 2 * length(u₀) + 1
    @test probes[1] == u₀
    h = max(sqrt(eps(eltype(u₀))), 1.0e-6)
    @test probes[2] ≈ u₀ + [h, 0]
    @test probes[3] ≈ u₀ - [h, 0]
    @test probes[4] ≈ u₀ + [0, h]
    @test probes[5] ≈ u₀ - [0, h]
end

@testset "_classify_noise_shape: CoupledSDEs" begin
    f_lin(u, p, t) = SA[-u[1], -u[2]]
    ds_iso = CoupledSDEs(f_lin, SA[0.0, 0.0]; noise_strength = 1.0)
    @test CriticalTransitions._classify_noise_shape(ds_iso) isa AdditiveNoise

    R = [cos(0.3) -sin(0.3); sin(0.3) cos(0.3)]
    Q_user = R * Diagonal([0.5, 2.0]) * R'
    ds_nd = CoupledSDEs(f_lin, SA[0.0, 0.0]; covariance = Q_user)
    @test CriticalTransitions._classify_noise_shape(ds_nd) isa GeneralNoise

    g1d(u, p, t) = SA[sqrt(1 + 0.3 * u[1]^2);;]
    f1d(u, p, t) = SA[-u[1]]
    ds_diagmult = CoupledSDEs(
        f1d, SA[1.0]; g = g1d, noise_prototype = SMatrix{1, 1}(0.0),
    )
    @test CriticalTransitions._classify_noise_shape(ds_diagmult) isa DiagonalNoise

    f2(u, p, t) = SA[-u[1], -u[2]]
    function g2(u, p, t)
        s11 = 1 + 0.2 * u[1]; s22 = 1 + 0.2 * u[2]
        s12 = 0.3 * u[2];     s21 = 0.3 * u[1]
        return @SMatrix [s11 s12; s21 s22]
    end
    ds_gen = CoupledSDEs(
        f2, SA[1.0, 0.0]; g = g2, noise_prototype = SMatrix{2, 2}(zeros(2, 2)),
    )
    @test CriticalTransitions._classify_noise_shape(ds_gen) isa GeneralNoise
end

@testset "_classify_noise_shape: rank-deficient rejection" begin
    function langevin(u, p, t)
        x, p_ = u
        return SA[p_, -x - 0.1 * p_]
    end
    g_langevin(u, p, t) = SA[0.0 0.0; 0.0 sqrt(0.2)]
    ds_langevin = CoupledSDEs(
        langevin, SA[0.0, 0.0]; g = g_langevin,
        noise_prototype = SMatrix{2, 2}(zeros(2, 2)),
    )
    err = try
        CriticalTransitions._classify_noise_shape(ds_langevin)
        nothing
    catch e
        e
    end
    @test err isa ArgumentError
    @test occursin("rank-deficient", err.msg)
end

@testset "proper_FW_system" begin
    f_lin(u, p, t) = SA[-u[1], -u[2]]
    ds_iso = CoupledSDEs(f_lin, SA[0.0, 0.0]; noise_strength = 1.0)
    @test CriticalTransitions.proper_FW_system(ds_iso) === nothing

    g1d(u, p, t) = SA[sqrt(1 + 0.3 * u[1]^2);;]
    f1d(u, p, t) = SA[-u[1]]
    ds_mult = CoupledSDEs(
        f1d, SA[1.0]; g = g1d, noise_prototype = SMatrix{1, 1}(0.0),
    )
    @test CriticalTransitions.proper_FW_system(ds_mult) === nothing
end
