using Test

@testset "auto-derived Hamiltonian derivative call forms" begin
    function linear_drift(u, p, t)
        return SA[-2 * u[1] + u[2], -u[1] - 3 * u[2]]
    end

    ds = CoupledSDEs(linear_drift, zeros(2); noise_strength = 0.5)
    sys = FreidlinWentzellHamiltonian(ds)

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
end
