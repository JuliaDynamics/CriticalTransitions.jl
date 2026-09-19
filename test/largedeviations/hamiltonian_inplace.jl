using Test

@testset "auto-derived Hamiltonian in-place evaluation" begin
    function linear_drift(u, p, t)
        return SA[-2 * u[1] + u[2], -u[1] - 3 * u[2]]
    end

    ds = CoupledSDEs(linear_drift, zeros(2); noise_strength = 0.5)
    sys = FreidlinWentzellHamiltonian(ds)
    @test sys isa FreidlinWentzellHamiltonian{true, 2}

    x = randn(2, 12)
    p = randn(2, 12)
    Hx_ref = sys.H_x(x, p)
    Hp_ref = sys.H_p(x, p)
    Hx = similar(x)
    Hp = similar(x)

    @test sys.H_x(Hx, x, p) === Hx
    @test sys.H_p(Hp, x, p) === Hp
    @test Hx ≈ Hx_ref
    @test Hp ≈ Hp_ref

    fill!(Hx, NaN)
    fill!(Hp, NaN)
    @test CT._eval_Hx!(Hx, sys, x, p) === Hx
    @test CT._eval_Hp!(Hp, sys, x, p) === Hp
    @test Hx ≈ Hx_ref
    @test Hp ≈ Hp_ref
end

@testset "buffered state-dependent diffusion evaluation" begin
    drift(u, p, t) = SA[-u[1], -u[2]]
    x = SA[0.2, -0.3]

    diagonal_noise(u, p, t) = SA[1 + 0.1 * u[1], 1 - 0.2 * u[2]]
    ds_diag = CoupledSDEs(drift, zeros(2); g = diagonal_noise)
    sys_diag = FreidlinWentzellHamiltonian(ds_diag)
    a_diag = zeros(2, 2)
    @test CT._eval_a!(a_diag, sys_diag.a, x) === a_diag
    @test a_diag ≈ Matrix(sys_diag.a(x))

    function coupled_noise(u, p, t)
        return @SMatrix [1 + 0.1 * u[1] 0.2 * u[2]; -0.15 * u[1] 1 - 0.1 * u[2]]
    end
    ds_coupled = CoupledSDEs(
        drift, zeros(2); g = coupled_noise,
        noise_prototype = SMatrix{2, 2}(zeros(2, 2)),
    )
    sys_coupled = FreidlinWentzellHamiltonian(ds_coupled)
    a_coupled = zeros(2, 2)
    @test CT._eval_a!(a_coupled, sys_coupled.a, x) === a_coupled
    @test a_coupled ≈ Matrix(sys_coupled.a(x))
end
