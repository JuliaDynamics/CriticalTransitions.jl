using CriticalTransitions
using Test
using LinearAlgebra

const CT = CriticalTransitions

function _underdamped_double_well(u, p, t)
    q, v = u
    return SA[v, q - q^3 - v]
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

@testset "MultipleShooting supports rank-deficient underdamped noise" begin
    ds = _rank_deficient_double_well()
    H = FreidlinWentzellHamiltonian(ds)
    init = _rank_deficient_initial_path()

    # GeometricGradient still uses an inverse diffusion metric and must reject this problem.
    @test_throws ArgumentError minimize_geometric_action(
        H, init, GeometricGradient(); maxiters = 1, show_progress = false
    )

    res = minimize_geometric_action(
        H,
        init,
        MultipleShooting(; nshoots = 10, maxiters = 200, abstol = 1.0e-8, reltol = 1.0e-7);
        show_progress = false,
    )

    # For U(q) = (q^2 - 1)^2/4 and fluctuation-dissipation noise in velocity,
    # the quasipotential barrier from (-1,0) to the saddle (0,0) is ΔU = 1/4.
    @test isapprox(res.action, 0.25; rtol = 3.0e-2)
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
        MultipleShooting(; nshoots = 10, maxiters = 200, abstol = 1.0e-8, reltol = 1.0e-7);
        show_progress = false,
    )

    @test isapprox(res.action, 0.25; rtol = 3.0e-2)
    @test _H_invariant_max(H, res) < 1.0e-5
end