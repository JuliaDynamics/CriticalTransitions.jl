using CriticalTransitions
using ModelingToolkit
using Test
using LinearAlgebra

@testset "ExtendedPhaseSpace KPO" begin
    λ = 3 / 1.21 * 2 / 295
    ω0 = 1.000
    ω = 1.000
    γ = 1 / 295
    η = 0
    α = -1

    function fu(u, v)
        return (-4 * γ * ω * u - 2 * λ * v - 4 * (ω0 - ω^2) * v - 3 * α * v * (u^2 + v^2)) /
               (8 * ω)
    end
    function fv(u, v)
        return (-4 * γ * ω * v - 2 * λ * u + 4 * (ω0 - ω^2) * u + 3 * α * u * (u^2 + v^2)) /
               (8 * ω)
    end

    dfvdv(u, v) = (-4 * γ * ω + 6 * α * u * v) / (8 * ω)
    dfudu(u, v) = (-4 * γ * ω - 6 * α * u * v) / (8 * ω)
    dfvdu(u, v) = (-2 * λ + 4 * (ω0 - ω^2) + 9 * α * u^2 + 3 * α * v^2) / (8 * ω)
    dfudv(u, v) = (-2 * λ - 4 * (ω0 - ω^2) - 3 * α * u^2 - 9 * α * v^2) / (8 * ω)

    function H_x(x, p) # ℜ² → ℜ²
        u, v = eachrow(x)
        pu, pv = eachrow(p)

        H_u = @. pu * dfudu(u, v) + pv * dfvdu(u, v)
        H_v = @. pu * dfudv(u, v) + pv * dfvdv(u, v)
        return Matrix([H_u H_v]')
    end
    function H_p(x, p) # ℜ² → ℜ²
        u, v = eachrow(x)
        pu, pv = eachrow(p)

        H_pu = @. pu + fu(u, v)
        H_pv = @. pv + fv(u, v)
        return Matrix([H_pu H_pv]')
    end

    @independent_variables t
    D = Differential(t)
    sts = @variables u(t) v(t)

    eqs = [D(u) ~ fu(u, v), D(v) ~ fv(u, v)]
    @mtkcompile sysMTK = System(eqs, t)
    prob = ODEProblem(sysMTK, Dict(sts .=> zeros(2)), (0.0, 100.0); jac=true)
    ds = CoupledODEs(prob)

    sys = ExtendedPhaseSpace{false,2}(H_x, H_p)
    sys′ = ExtendedPhaseSpace(ds)

    Nt = 500  # number of discrete time steps
    p_r = rand(2, Nt)
    x_r = rand(2, Nt)

    # MTK v10+ reorders variables: unknowns are [v, u] instead of [u, v]
    # Swap rows to match MTK's internal ordering
    x_r_swapped = [x_r[2, :]'; x_r[1, :]']
    p_r_swapped = [p_r[2, :]'; p_r[1, :]']

    result_H_x = sys′.H_x(x_r_swapped, p_r_swapped)
    result_H_p = sys′.H_p(x_r_swapped, p_r_swapped)

    # Swap results back to compare with manual system
    @test [result_H_x[2, :]'; result_H_x[1, :]'] ≈ sys.H_x(x_r, p_r)
    @test [result_H_p[2, :]'; result_H_p[1, :]'] ≈ sys.H_p(x_r, p_r)
end

@testset "ExtendedPhaseSpace MTK" begin
    @independent_variables t
    D = Differential(t)
    sts = @variables u(t) v(t)

    @parameters λ = 3 / 1.21 * 2 / 295 ω0 = 1.0 ω = 1.0 γ = 1 / 295 η = 0 α = -1

    eqs = [
        D(u) ~
        (-4 * γ * ω * u - 2 * λ * v - 4 * (ω0 - ω^2) * v - 3 * α * v * (u^2 + v^2)) /
        (8 * ω),
        D(v) ~
        (-4 * γ * ω * v - 2 * λ * u + 4 * (ω0 - ω^2) * u + 3 * α * u * (u^2 + v^2)) /
        (8 * ω),
    ]
    @mtkcompile sysMTK = System(eqs, t)
    prob = ODEProblem(sysMTK, Dict(sts .=> zeros(2)), (0.0, 100.0); jac=true)
    ds = CoupledODEs(prob)
    sys = ExtendedPhaseSpace(ds)

    @test sys.H_x(zeros(2), zeros(2)) ≈ zeros(2)
    @test sys.H_p(zeros(2), zeros(2)) ≈ zeros(2)
end

@testset "sgMAM GeometricGradient" begin
    H_x(x, p) = zeros(size(x))
    H_p(x, p) = ones(size(x))
    sys = ExtendedPhaseSpace{false,2}(H_x, H_p)

    xx = collect(range(-1.0, 1.0; length=20))
    yy = 0.3 .* (-xx .^ 2 .+ 1)
    x_initial = Matrix([xx yy]')

    res_small = minimize_simple_geometric_action(
        sys, x_initial, GeometricGradient(; stepsize=1e-6); maxiters=2, show_progress=false
    )
    res_large = minimize_simple_geometric_action(
        sys, x_initial, GeometricGradient(; stepsize=1.0); maxiters=2, show_progress=false
    )

    @test res_small.action != res_large.action
end
