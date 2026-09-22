using CriticalTransitions, StaticArrays
using CriticalTransitions.CTLibrary: ou_multiplicative_1d, linear_offdiag_2d_sde
using Test
using LinearAlgebra
using Random

const _make_1d_ou = ou_multiplicative_1d
const _make_2d_offdiag = linear_offdiag_2d_sde

@testset "sgMAM end-to-end: 1D OU multiplicative converges" begin
    sys = FreidlinWentzellHamiltonian(_make_1d_ou(0.3))
    Nt = 80
    x_initial = reshape(collect(range(1.0, -1.0; length = Nt)), 1, Nt)
    res = minimize_geometric_action(
        sys, x_initial, GeometricGradient(; stepsize = 1.0);
        maxiters = 500, show_progress = false,
    )
    @test isfinite(res.action)
end

@testset "sgMAM end-to-end: 2D off-diagonal multiplicative converges" begin
    Random.seed!(0)
    sys = FreidlinWentzellHamiltonian(_make_2d_offdiag())
    Nt = 60
    xx = collect(range(1.0, 0.0; length = Nt))
    yy = collect(range(0.0, 1.0; length = Nt))
    x_initial = Matrix([xx yy]')
    res = minimize_geometric_action(
        sys, x_initial, GeometricGradient(; stepsize = 1.0);
        maxiters = 300, show_progress = false,
    )
    @test isfinite(res.action)
end

@testset "fw_action: 1D OU multiplicative vs analytic Simpson" begin
    α = 0.3
    ds = _make_1d_ou(α)

    N, T = 200, 1.0
    path = reduce(hcat, range([1.0], [0.0]; length = N))
    time = range(0.0, T; length = N)
    S = fw_action(ds, path, time)

    function simpson(f, a, b, n)
        h = (b - a) / n; s = f(a) + f(b)
        for i in 1:2:(n - 1)
            s += 4 * f(a + i * h)
        end
        for i in 2:2:(n - 2)
            s += 2 * f(a + i * h)
        end
        return s * h / 3
    end
    s_norm = 1 + α
    integrand = t -> t^2 / (1 + α * (1 - t)^2)
    analytic = s_norm * simpson(integrand, 0, 1, 10_000) / 2
    @test isapprox(S, analytic; rtol = 1.0e-4)
end

@testset "gMAM diagonal multiplicative converges" begin
    Random.seed!(0)
    ds = _make_1d_ou(0.3)
    Nt = 80
    x_initial = reshape(collect(range(1.0, -1.0; length = Nt)), 1, Nt)
    res_g = minimize_geometric_action(
        ds, x_initial, GeometricGradient(; stepsize = 1.0);
        maxiters = 500, show_progress = false,
    )
    @test isfinite(res_g.action)
end

@testset "MAM (FW) on multiplicative noise agrees with gMAM at large T" begin
    # Cross-algorithm correctness check: MAM (time-parameterized FW) and gMAM
    # (geometric / time-eliminated FW) compute the same physical instanton
    # action up to discretization. The continuous identity
    # `inf_T S_T^{FW}[φ] = S_geo[φ]` is exact, but the discrete trapezoid is
    # only reparameterization-invariant up to `O(1/N²)` curvature terms, so
    # we expect a few-percent gap at moderate `N`, not bitwise agreement.
    Random.seed!(0)
    ds = _make_1d_ou(0.3)
    Nt = 80

    xinit_g = reshape(collect(range(1.0, -1.0; length = Nt)), 1, Nt)
    res_g = minimize_geometric_action(
        ds, xinit_g, GeometricGradient(; stepsize = 1.0);
        maxiters = 2000, show_progress = false,
    )
    S_g = res_g.action

    # Run MAM (FW) over a range of T. At small T the uniform-Δt constraint is
    # restrictive (action is dominated by the kinetic `|φ̇|²/T` term); as T
    # grows the optimizer can cluster path points to simulate the
    # FW-natural non-uniform speed profile, and the action drops toward the
    # gMAM value.
    init = reduce(hcat, range([1.0], [-1.0]; length = Nt))
    Ts = [0.5, 1.0, 4.0, 16.0]
    S_mam = map(Ts) do T
        res = minimize_action(ds, init, T; maxiters = 2000, show_progress = false)
        return res.action
    end

    @test issorted(S_mam; rev = true)
    @test S_mam[1] > 3 * S_g
    @test S_mam[end] ≈ S_g rtol = 0.05

    # `functional = "OM"` is rejected for multiplicative noise (the OM
    # correction term is only implemented for additive diffusion).
    @test_throws ArgumentError minimize_action(
        ds, init, 1.0;
        functional = "OM", noise_strength = 0.1,
        maxiters = 5, show_progress = false,
    )
end

@testset "sgMAM additive non-diagonal covariance converges" begin
    function meier_stein(u, p, t)
        x, y = u
        return SA[x - x^3 - 10 * x * y^2, -(1 + x^2) * y]
    end
    Nt = 60
    xx = range(-1.0, 1.0; length = Nt)
    yy = 0.3 .* (-xx .^ 2 .+ 1)
    x_initial = Matrix([xx yy]')

    θ = 0.4
    R = [cos(θ) -sin(θ); sin(θ) cos(θ)]
    D = Diagonal([0.5, 2.0])
    Q_rot = R * D * R'
    ds_rot = CoupledSDEs(meier_stein, zeros(2); covariance = Q_rot)
    sys_rot = FreidlinWentzellHamiltonian(ds_rot)

    res = minimize_geometric_action(
        sys_rot, x_initial, GeometricGradient(; stepsize = 1.0);
        maxiters = 500, show_progress = false,
    )
    @test isfinite(res.action)
end

@testset "gMAM general multiplicative converges" begin
    Random.seed!(0)
    ds = _make_2d_offdiag()
    Nt = 60
    xx = collect(range(1.0, 0.0; length = Nt))
    yy = collect(range(0.0, 1.0; length = Nt))
    x_initial = Matrix([xx yy]')
    res_g = minimize_geometric_action(
        ds, x_initial, GeometricGradient(; stepsize = 1.0);
        maxiters = 300, show_progress = false,
    )
    @test isfinite(res_g.action)
end

# Regression: the state-dependent gMAM/sgMAM branches must reduce to the additive
# branch when `a(x) ≡ I`. We construct a state-dependent diffusion `σ(x) = R(x)`
# whose product `σσᵀ = R Rᵀ ≡ I`, so classification routes it through the
# DiagonalNoise / GeneralNoise paths while the resulting Hamiltonian is identical
# to the constant-`a` case. The converged action should match an additive-`σ ≡ I`
# system on the same drift.
@testset "gMAM state-dep a ≡ I matches additive a ≡ I" begin
    Random.seed!(0)
    function ms_drift(u, p, t)
        x, y = u
        return SA[x - x^3 - 10 * x * y^2, -(1 + x^2) * y]
    end

    # Reference: additive identity noise.
    ds_add = CoupledSDEs(ms_drift, zeros(2); noise_strength = 1.0)

    # State-dependent σ(x) = R(θ(x)) with θ(x) = 0.5 * x[1]; σσᵀ = R Rᵀ = I.
    function g_rotation(u, p, t)
        c = cos(0.5 * u[1]); s = sin(0.5 * u[1])
        return @SMatrix [c -s; s c]
    end
    ds_sd = CoupledSDEs(
        ms_drift, zeros(2); g = g_rotation,
        noise_prototype = SMatrix{2, 2}(zeros(2, 2)),
    )

    Nt = 40
    xx = range(-1.0, 1.0; length = Nt)
    yy = 0.3 .* (-xx .^ 2 .+ 1)
    x_initial = Matrix([xx yy]')

    res_add = minimize_geometric_action(
        ds_add, x_initial, GeometricGradient(; stepsize = 1.0);
        maxiters = 500, show_progress = false,
    )
    res_sd = minimize_geometric_action(
        ds_sd, x_initial, GeometricGradient(; stepsize = 1.0);
        maxiters = 500, show_progress = false,
    )

    @test isapprox(res_add.action, res_sd.action; rtol = 1.0e-6)
end
