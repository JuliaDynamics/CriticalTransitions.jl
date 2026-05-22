using CriticalTransitions, StaticArrays
using Test
using LinearAlgebra

const CT = CriticalTransitions

@testset "minimize_geometric_action dispatches: CoupledSDEs (gMAM) ≈ FreidlinWentzellHamiltonian (sgMAM)" begin
    function meier_stein(u, p, t)
        x, y = u
        return SA[x - x^3 - 10 * x * y^2, -(1 + x^2) * y]
    end
    σ = 0.25
    ds = CoupledSDEs(meier_stein, zeros(2); noise_strength = σ)
    sys_h = FreidlinWentzellHamiltonian(ds)

    Nt = 60
    xx = range(-1.0, 1.0; length = Nt)
    yy = 0.3 .* (-xx .^ 2 .+ 1)
    x_initial = Matrix([xx yy]')

    res_g = minimize_geometric_action(
        ds, x_initial, GeometricGradient(; stepsize = 1.0);
        maxiters = 500, show_progress = false,
    )
    res_s = minimize_geometric_action(
        sys_h, x_initial, GeometricGradient(; stepsize = 1.0);
        maxiters = 500, show_progress = false,
    )
    @test isapprox(res_g.action, res_s.action; rtol = 1.0e-2)
end

@testset "minimize_geometric_action rejects CoupledODEs (#263)" begin
    f_lin(u, p, t) = SA[-u[1], -u[2]]
    ds_ode = CoupledODEs(f_lin, SA[0.0, 0.0])
    xx = collect(range(-1.0, 1.0; length = 20))
    yy = 0.3 .* (-xx .^ 2 .+ 1)
    path = Matrix([xx yy]')
    @test_throws MethodError minimize_geometric_action(ds_ode, path)
end

@testset "minimize_simple_geometric_action is removed (#326)" begin
    @test !isdefined(CriticalTransitions, :minimize_simple_geometric_action)
end
