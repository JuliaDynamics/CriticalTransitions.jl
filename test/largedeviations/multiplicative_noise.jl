using CriticalTransitions, StaticArrays
using Test
using LinearAlgebra
using Random

const CT = CriticalTransitions

# Simple Simpson's rule for analytic baselines, avoiding an extra test dep
function _simpson(f, a, b, n)
    h = (b - a) / n
    s = f(a) + f(b)
    for i in 1:2:(n - 1)
        s += 4 * f(a + i * h)
    end
    for i in 2:2:(n - 2)
        s += 2 * f(a + i * h)
    end
    return s * h / 3
end

@testset "fw_action: state-dependent 1D OU vs analytic" begin
    # SDE: dX = -X dt + sqrt(1 + α X²) dW;  a(x) = 1 + α x²
    # Straight path x(t) = 1 - t over T = 1
    α = 0.3
    b1d(u, p, t) = SA[-u[1]]
    g1d(u, p, t) = SA[sqrt(1 + α * u[1]^2);;]
    ds = CoupledSDEs(b1d, SA[1.0]; g = g1d, noise_prototype = SMatrix{1, 1}(0.0))

    N, T = 200, 1.0
    path = reduce(hcat, range([1.0], [0.0]; length = N))
    time = range(0.0, T; length = N)
    S = fw_action(ds, path, time)

    # Convention B: the FW rate function uses a_shape(x) = a(x)/s with s = a(u₀) at the
    # system's reference state u₀ = 1.0. Here a(x) = 1 + α x², so s = 1 + α. Analytic:
    #   S = (1/2) ∫₀¹ t² · s / (1 + α (1-t)²) dt  =  s · (1/2) ∫₀¹ t² / (1+α(1-t)²) dt.
    integrand = t -> t^2 / (1 + α * (1 - t)^2)
    s = 1 + α
    analytic = s * _simpson(integrand, 0, 1, 10_000) / 2

    @test isapprox(S, analytic; rtol = 1.0e-4)
end

@testset "geometric_action: additive limit recovers constant-Q value" begin
    # When σ(x) is constant, the auto-dispatched multiplicative path must agree
    # with the additive path (after matching the additive-side normalization).
    p = [0.1, 3, 1, 1, 1, 0]
    σ = 0.1
    ds_add = CoupledSDEs(CT.CTLibrary.fitzhugh_nagumo, zeros(2), p; noise_strength = σ)

    x_i = SA[sqrt(2 / 3), sqrt(2 / 27)]
    x_f = SA[0.001, 0.0]
    N = 60
    path = reduce(hcat, range(x_i, x_f; length = N))

    S_add = geometric_action(ds_add, path)

    # Match the additive normalization in a manually-built multiplicative system.
    Q = covariance_matrix(ds_add)
    Q_norm = CT.normalize_covariance!(Matrix(Q))
    σ_const = sqrt(Q_norm)
    g_const(u, p, t) = SMatrix{2, 2}(σ_const)
    ds_mul = CoupledSDEs(
        CT.CTLibrary.fitzhugh_nagumo, zeros(2), p;
        g = g_const, noise_prototype = SMatrix{2, 2}(zeros(2, 2)),
    )

    S_mul = geometric_action(ds_mul, path)
    @test isapprox(S_add, S_mul; rtol = 1.0e-10)
end

@testset "minimize_action via Optimization.jl (multiplicative auto-dispatch)" begin
    # With a constant σ that already matches the additive-side normalization,
    # the multiplicative path must coincide with the additive minimizer.
    Random.seed!(0)
    σ = 0.18
    ou_add = CoupledSDEs((u, p, t) -> -u, SA[1.0]; noise_strength = σ)
    x0 = -1.0
    xT = 2.0
    T = 10.0
    N = 51
    res_ref = minimize_action(
        ou_add, SA[x0], SA[xT], T; npoints = N, maxiters = 2000, show_progress = false
    )

    Q = covariance_matrix(ou_add)
    Q_norm = CT.normalize_covariance!(Matrix(Q))
    σ_const = sqrt(Q_norm[1, 1])
    g_const(u, p, t) = SA[σ_const;;]
    ou_mul = CoupledSDEs(
        (u, p, t) -> -u, SA[1.0]; g = g_const, noise_prototype = SMatrix{1, 1}(0.0)
    )

    times = range(0.0, T; length = N)
    S_with_a = fw_action(ou_mul, Matrix(res_ref.path)', times)
    S_ref = fw_action(ou_add, Matrix(res_ref.path)', times)

    @test isapprox(S_with_a, S_ref; rtol = 1.0e-10)
end

@testset "gMAM (GeometricGradient) diagonal multiplicative vs Adam" begin
    using OptimizationOptimisers
    Random.seed!(0)
    α = 0.3
    b1d(u, p, t) = SA[-u[1]]
    g1d(u, p, t) = SA[sqrt(1 + α * u[1]^2);;]
    ds = CoupledSDEs(b1d, SA[1.0]; g = g1d, noise_prototype = SMatrix{1, 1}(0.0))

    N = 80
    x_initial = reshape(collect(range(1.0, -1.0; length = N)), 1, N)

    res_g = minimize_geometric_action(
        ds, x_initial, GeometricGradient(; stepsize = 1.0);
        maxiters = 500, show_progress = false,
    )
    res_ref = minimize_geometric_action(
        ds, x_initial, OptimizationOptimisers.Adam(0.01);
        maxiters = 3000, show_progress = false,
    )
    @test isapprox(res_g.action, res_ref.action; rtol = 0.05)
end

@testset "gMAM (GeometricGradient) general multiplicative vs Adam" begin
    using OptimizationOptimisers
    Random.seed!(0)
    b2(u, p, t) = SA[-u[1], -u[2]]
    function g2(u, p, t)
        s11 = 1 + 0.2 * u[1]
        s22 = 1 + 0.2 * u[2]
        s12 = 0.3 * u[2]
        s21 = 0.3 * u[1]
        return @SMatrix [s11 s12; s21 s22]
    end
    ds = CoupledSDEs(b2, SA[1.0, 0.0]; g = g2, noise_prototype = SMatrix{2, 2}(zeros(2, 2)))

    N = 60
    xx = collect(range(1.0, 0.0; length = N))
    yy = collect(range(0.0, 1.0; length = N))
    x_initial = Matrix([xx yy]')

    res_g = minimize_geometric_action(
        ds, x_initial, GeometricGradient(; stepsize = 1.0);
        maxiters = 300, show_progress = false,
    )
    res_ref = minimize_geometric_action(
        ds, x_initial, OptimizationOptimisers.Adam(0.01);
        maxiters = 3000, show_progress = false,
    )
    @test isfinite(res_g.action)
    @test isapprox(res_g.action, res_ref.action; rtol = 0.1)
end

@testset "sgMAM diagonal multiplicative vs Adam reference" begin
    using OptimizationOptimisers
    Random.seed!(0)
    α = 0.3
    b1d(u, p, t) = SA[-u[1]]
    g1d(u, p, t) = SA[sqrt(1 + α * u[1]^2);;]
    ds = CoupledSDEs(b1d, SA[1.0]; g = g1d, noise_prototype = SMatrix{1, 1}(0.0))

    sys = FreidlinWentzellHamiltonian(ds)
    @test sys.a_func !== nothing
    @test sys.a_func([0.5]) isa AbstractVector
    # `a_func` returns the *shape* a(x)/s with s = a(u₀) at the reference state u₀ = 1.0.
    # At x = 0.5: a(0.5)/s = (1 + α · 0.25) / (1 + α).
    @test sys.a_func([0.5])[1] ≈ (1 + α * 0.5^2) / (1 + α)

    N = 80
    x_initial = reshape(collect(range(1.0, -1.0; length = N)), 1, N)

    res_sg = minimize_simple_geometric_action(
        sys, x_initial, GeometricGradient(; stepsize = 1.0);
        maxiters = 500, show_progress = false,
    )
    res_ref = minimize_geometric_action(
        ds, x_initial, OptimizationOptimisers.Adam(0.01);
        maxiters = 3000, show_progress = false,
    )
    @test isapprox(res_sg.action, res_ref.action; rtol = 0.05)
end

@testset "sgMAM general (non-diagonal) multiplicative noise" begin
    using OptimizationOptimisers
    Random.seed!(0)
    b2(u, p, t) = SA[-u[1], -u[2]]
    function g2(u, p, t)
        s11 = 1 + 0.2 * u[1]
        s22 = 1 + 0.2 * u[2]
        s12 = 0.3 * u[2]
        s21 = 0.3 * u[1]
        return @SMatrix [s11 s12; s21 s22]
    end
    ds = CoupledSDEs(b2, SA[1.0, 0.0]; g = g2, noise_prototype = SMatrix{2, 2}(zeros(2, 2)))

    sys = FreidlinWentzellHamiltonian(ds)
    @test sys.a_func([0.5, 0.3]) isa AbstractMatrix
    @test !LinearAlgebra.isdiag(sys.a_func([0.5, 0.3]))

    N = 60
    xx = collect(range(1.0, 0.0; length = N))
    yy = collect(range(0.0, 1.0; length = N))
    x_initial = Matrix([xx yy]')

    res_sg = minimize_simple_geometric_action(
        sys, x_initial, GeometricGradient(; stepsize = 1.0);
        maxiters = 300, show_progress = false,
    )
    @test isfinite(res_sg.action)

    res_ref = minimize_geometric_action(
        ds, x_initial, OptimizationOptimisers.Adam(0.01);
        maxiters = 3000, show_progress = false,
    )
    @test isapprox(res_sg.action, res_ref.action; rtol = 0.1)
end

@testset "Multiplicative → additive limit (1D OU, sgMAM and gMAM)" begin
    # For dX = -X dt + dW the geometric FW action from (1) to (-1) is
    # 2·V(-1) = 1 (downhill 1 → 0 is free, uphill 0 → -1 costs 1). As α → 0,
    # the multiplicative system dX = -X dt + √(1+αX²) dW must recover the same value.
    additive(u, p, t) = SA[-u[1]]
    ds_add = CoupledSDEs(additive, SA[1.0]; noise_strength = 1.0)

    α_small = 1.0e-4
    b1d(u, p, t) = SA[-u[1]]
    g1d(u, p, t) = SA[sqrt(1 + α_small * u[1]^2);;]
    ds_mul = CoupledSDEs(
        b1d, SA[1.0]; g = g1d, noise_prototype = SMatrix{1, 1}(0.0)
    )

    N = 100
    x_initial = reshape(collect(range(1.0, -1.0; length = N)), 1, N)

    # gMAM: additive analytic value, multiplicative converges to additive as α → 0.
    res_add_g = minimize_geometric_action(
        ds_add, x_initial, GeometricGradient(; stepsize = 1.0);
        maxiters = 500, show_progress = false,
    )
    res_mul_g = minimize_geometric_action(
        ds_mul, x_initial, GeometricGradient(; stepsize = 1.0);
        maxiters = 500, show_progress = false,
    )
    @test isapprox(res_add_g.action, 1.0; atol = 3.0e-2)
    @test isapprox(res_mul_g.action, res_add_g.action; atol = 5.0e-3)

    # sgMAM: same check via the Hamiltonian path.
    sys_add = FreidlinWentzellHamiltonian(ds_add)
    sys_mul = FreidlinWentzellHamiltonian(ds_mul)
    res_add_sg = minimize_simple_geometric_action(
        sys_add, x_initial, GeometricGradient(; stepsize = 1.0);
        maxiters = 500, show_progress = false,
    )
    res_mul_sg = minimize_simple_geometric_action(
        sys_mul, x_initial, GeometricGradient(; stepsize = 1.0);
        maxiters = 500, show_progress = false,
    )
    @test isapprox(res_add_sg.action, 1.0; atol = 3.0e-2)
    @test isapprox(res_mul_sg.action, res_add_sg.action; atol = 5.0e-3)
end
