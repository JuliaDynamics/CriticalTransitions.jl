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

_uphill_rhs(y) = SA[y[2], y[1] - y[1]^3 + y[2]]

function _rk4_uphill_step(y, h)
    k1 = _uphill_rhs(y)
    k2 = _uphill_rhs(y + (h / 2) * k1)
    k3 = _uphill_rhs(y + (h / 2) * k2)
    k4 = _uphill_rhs(y + h * k3)
    return y + (h / 6) * (k1 + 2k2 + 2k3 + k4)
end

# Exact equilibrium instanton in configuration space: integrate the anti-damped uphill
# dynamics backward from the saddle's stable direction, reverse the heteroclinic, then
# resample it uniformly in Euclidean arclength. This is a test oracle / admissible warm
# start, not a production special case.
function _rank_deficient_admissible_path(N = 100)
    nsteps = 5000
    h = -0.01
    δ = 1.0e-6
    λstable = (1 - sqrt(5.0)) / 2
    y = SA[-δ, -λstable * δ]

    backward = Matrix{Float64}(undef, 2, nsteps + 1)
    backward[:, 1] .= y
    for i in 1:nsteps
        y = _rk4_uphill_step(y, h)
        backward[:, i + 1] .= y
    end

    orbit = Matrix{Float64}(undef, 2, nsteps + 3)
    orbit[:, 1] .= (-1.0, 0.0)
    orbit[:, 2:(end - 1)] .= backward[:, end:-1:1]
    orbit[:, end] .= (0.0, 0.0)

    arclength = zeros(Float64, size(orbit, 2))
    for i in 2:length(arclength)
        arclength[i] = arclength[i - 1] + norm(orbit[:, i] - orbit[:, i - 1])
    end

    out = Matrix{Float64}(undef, 2, N)
    targets = range(0.0, arclength[end]; length = N)
    j = 1
    for (i, target) in enumerate(targets)
        while j < length(arclength) - 1 && arclength[j + 1] < target
            j += 1
        end
        if target == arclength[end]
            out[:, i] .= orbit[:, end]
        else
            Δs = arclength[j + 1] - arclength[j]
            θ = Δs > 0 ? (target - arclength[j]) / Δs : 0.0
            out[:, i] .= (1 - θ) .* orbit[:, j] .+ θ .* orbit[:, j + 1]
        end
    end
    out[:, 1] .= (-1.0, 0.0)
    out[:, end] .= (0.0, 0.0)
    return out
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

@testset "Rank-deficient underdamped Hamiltonian oracle" begin
    H = FreidlinWentzellHamiltonian(_rank_deficient_double_well())
    for q in (-0.9, -0.6, -0.3), v in (-0.2, 0.15)
        x = [q, v]
        p = [q^3 - q, v]
        @test abs(CT._hamiltonian_value(H, x, p)) < 1.0e-12
        Hp = H.H_p(reshape(x, 2, 1), reshape(p, 2, 1))[:, 1]
        @test Hp ≈ [v, q - q^3 + v] atol = 1.0e-12
    end
end

@testset "MultipleShooting supports rank-deficient underdamped noise" begin
    ds = _rank_deficient_double_well()
    H = FreidlinWentzellHamiltonian(ds)
    straight = _rank_deficient_initial_path()
    init = _rank_deficient_admissible_path()

    # GeometricGradient still uses an inverse diffusion metric and must reject this problem.
    @test_throws ArgumentError minimize_geometric_action(
        H, straight, GeometricGradient(); maxiters = 1, show_progress = false
    )

    res = minimize_geometric_action(
        H,
        init,
        MultipleShooting(
            ; nshoots = 10, nlsolve = _rank_deficient_trust_region(), maxiters = 200,
            abstol = 1.0e-8, reltol = 1.0e-7,
        );
        show_progress = false,
    )

    # The package trace-normalizes raw rank-one velocity noise to a = diag(0, 2), for which
    # V(q,v) = U(q) + v²/2 and the barrier from (-1,0) to (0,0) is ΔU = 1/4.
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
    init = R * _rank_deficient_admissible_path()

    res = minimize_geometric_action(
        H,
        init,
        MultipleShooting(
            ; nshoots = 10, nlsolve = _rank_deficient_trust_region(), maxiters = 200,
            abstol = 1.0e-8, reltol = 1.0e-7,
        );
        show_progress = false,
    )

    @test isapprox(res.action, 0.25; rtol = 3.0e-2)
    @test _H_invariant_max(H, res) < 1.0e-5
end
