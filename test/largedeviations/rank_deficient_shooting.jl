using CriticalTransitions
using Test
using LinearAlgebra

const CT = CriticalTransitions
const Γ_RD = 3.0

function _underdamped_double_well(u, p, t)
    q, v = u
    return SA[v, q - q^3 - Γ_RD * v]
end

_rank_one_velocity_noise(u, p, t) = SA[0.0 0.0; 0.0 1.0]

function _rank_deficient_double_well()
    return CoupledSDEs(
        _underdamped_double_well,
        SA[-1.0, 0.0];
        g = _rank_one_velocity_noise,
        noise_prototype = SMatrix{2, 2}(zeros(2, 2)),
    )
end

function _rank_deficient_initial_path(N = 30)
    q = collect(range(-1.0, 0.0; length = N))
    return Matrix(vcat(q', zeros(1, N)))
end

function _H_invariant_max(H, res)
    D = size(res.generalized_momentum, 2)
    return maximum(
        abs(
            CT._hamiltonian_value(
                H,
                [res.path[i][k] for k in 1:D],
                [res.generalized_momentum[i, k] for k in 1:D],
            ),
        ) for i in eachindex(res.path)
    )
end

function _rank_deficient_trust_region()
    return CT.NonlinearSolveFirstOrder.TrustRegion(
        ; autodiff = CT.AutoForwardDiff(),
        linsolve = CT.LinearSolve.UMFPACKFactorization(),
    )
end

@testset "Rank-deficient second-order Langevin Hamiltonian oracle" begin
    H = FreidlinWentzellHamiltonian(_rank_deficient_double_well())
    for q in (-0.9, -0.6, -0.3), v in (-0.2, 0.15)
        x = [q, v]
        p = Γ_RD .* [q^3 - q, v]
        @test abs(CT._hamiltonian_value(H, x, p)) < 1.0e-12
        Hp = H.H_p(reshape(x, 2, 1), reshape(p, 2, 1))[:, 1]
        @test Hp ≈ [v, q - q^3 + Γ_RD * v] atol = 1.0e-12
    end
end

@testset "MultipleShooting supports rank-deficient second-order Langevin noise" begin
    ds = _rank_deficient_double_well()
    H = FreidlinWentzellHamiltonian(ds)
    init = _rank_deficient_initial_path()

    # GeometricGradient still uses an inverse diffusion metric and must reject this problem.
    @test_throws ArgumentError minimize_geometric_action(
        H, init, GeometricGradient(); maxiters = 1, show_progress = false
    )

    # Diagnostic coarse solve: if two segments converge, singular shooting should be
    # initialized by segment-count continuation rather than a one-shot 10-segment solve.
    res = minimize_geometric_action(
        H,
        init,
        MultipleShooting(
            ; nshoots = 2, nlsolve = _rank_deficient_trust_region(), maxiters = 200,
            abstol = 1.0e-8, reltol = 1.0e-7,
        );
        show_progress = false,
    )

    # Trace normalization gives a = diag(0, 2). For damping Γ, the exact quasipotential is
    # V(q,v) = Γ [U(q) + v²/2], hence S((-1,0) → (0,0)) = Γ ΔU = 3/4.
    @test isapprox(res.action, 0.75; rtol = 3.0e-2)
    @test _H_invariant_max(H, res) < 1.0e-5
    @test isapprox(res.path[1][1], -1.0; atol = 1.0e-5)
    @test isapprox(res.path[end][1], 0.0; atol = 1.0e-4)
end

@testset "Rank-deficient shooting is rotation invariant" begin
    θ = π / 5
    c, s = cos(θ), sin(θ)
    R = SA[c -s; s c]
    RT = R'

    function drift_rot(u, p, t)
        y = RT * u
        return R * _underdamped_double_well(y, p, t)
    end

    σ0 = SA[0.0 0.0; 0.0 1.0]
    σrot = R * σ0
    noise_rot(u, p, t) = σrot

    xa = R * SA[-1.0, 0.0]
    ds = CoupledSDEs(
        drift_rot,
        xa;
        g = noise_rot,
        noise_prototype = SMatrix{2, 2}(zeros(2, 2)),
    )
    H = FreidlinWentzellHamiltonian(ds)
    init = R * _rank_deficient_initial_path()

    res = minimize_geometric_action(
        H,
        init,
        MultipleShooting(
            ; nshoots = 2, nlsolve = _rank_deficient_trust_region(), maxiters = 200,
            abstol = 1.0e-8, reltol = 1.0e-7,
        );
        show_progress = false,
    )

    @test isapprox(res.action, 0.75; rtol = 3.0e-2)
    @test _H_invariant_max(H, res) < 1.0e-5
end
