using CriticalTransitions
using Test
using LinearAlgebra


const CT_ns = CriticalTransitions

@testset "_check_rank!" begin
    CT_ns._check_rank!(Matrix{Float64}(LinearAlgebra.I(3)))
    CT_ns._check_rank!(LinearAlgebra.Diagonal([1.0, 2.0, 3.0]))
    @test_throws ArgumentError CT_ns._check_rank!(LinearAlgebra.Diagonal([1.0, 0.0, 3.0]))
    @test_throws ArgumentError CT_ns._check_rank!([1.0 1.0; 1.0 1.0])
end

@testset "_isdiag_numerical" begin
    @test CT_ns._isdiag_numerical(LinearAlgebra.Diagonal([1.0, 2.0]))
    @test CT_ns._isdiag_numerical([1.0 0.0; 0.0 2.0])
    @test !CT_ns._isdiag_numerical([1.0 0.3; 0.3 2.0])
end

@testset "_validate_and_classify_a: CoupledSDEs" begin
    f_lin(u, p, t) = SA[-u[1], -u[2]]

    ds_iso = CoupledSDEs(f_lin, SA[0.0, 0.0]; noise_strength = 1.0)
    a_iso = CT_ns._trace_normalized_a(ds_iso)
    @test CT_ns._validate_and_classify_a(a_iso, [0.0, 0.0]) === true

    R = [cos(0.3) -sin(0.3); sin(0.3) cos(0.3)]
    Q_user = R * Diagonal([0.5, 2.0]) * R'
    ds_nd = CoupledSDEs(f_lin, SA[0.0, 0.0]; covariance = Q_user)
    a_nd = CT_ns._trace_normalized_a(ds_nd)
    @test CT_ns._validate_and_classify_a(a_nd, [0.0, 0.0]) === false

    g1d(u, p, t) = SA[sqrt(1 + 0.3 * u[1]^2);;]
    f1d(u, p, t) = SA[-u[1]]
    ds_diagmult = CoupledSDEs(
        f1d, SA[1.0]; g = g1d, noise_prototype = SMatrix{1, 1}(0.0),
    )
    a_dm = CT_ns._trace_normalized_a(ds_diagmult)
    @test CT_ns._validate_and_classify_a(a_dm, [1.0]) === true

    f2(u, p, t) = SA[-u[1], -u[2]]
    function g2(u, p, t)
        s11 = 1 + 0.2 * u[1]; s22 = 1 + 0.2 * u[2]
        s12 = 0.3 * u[2];     s21 = 0.3 * u[1]
        return @SMatrix [s11 s12; s21 s22]
    end
    ds_gen = CoupledSDEs(
        f2, SA[1.0, 0.0]; g = g2, noise_prototype = SMatrix{2, 2}(zeros(2, 2)),
    )
    a_g = CT_ns._trace_normalized_a(ds_gen)
    @test CT_ns._validate_and_classify_a(a_g, [1.0, 0.0]) === false
end

@testset "_validate_and_classify_a: rank-deficient rejection" begin
    function langevin(u, p, t)
        x, p_ = u
        return SA[p_, -x - 0.1 * p_]
    end
    g_langevin(u, p, t) = SA[0.0 0.0; 0.0 sqrt(0.2)]
    ds_langevin = CoupledSDEs(
        langevin, SA[0.0, 0.0]; g = g_langevin,
        noise_prototype = SMatrix{2, 2}(zeros(2, 2)),
    )
    a = CT_ns._trace_normalized_a(ds_langevin)
    err = try
        CT_ns._validate_and_classify_a(a, [0.0, 0.0])
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
    @test CT_ns.proper_FW_system(ds_iso) === nothing

    g1d(u, p, t) = SA[sqrt(1 + 0.3 * u[1]^2);;]
    f1d(u, p, t) = SA[-u[1]]
    ds_mult = CoupledSDEs(
        f1d, SA[1.0]; g = g1d, noise_prototype = SMatrix{1, 1}(0.0),
    )
    @test CT_ns.proper_FW_system(ds_mult) === nothing
end
