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

# Diagnostic only: obtain a globally admissible phase-space seed by solving nearby
# full-rank Hamiltonians and continuing the regularization toward the singular problem.
# The physical shooting equations and final action are always evaluated at δ = 0.
function _regularized_hamiltonian(H, xref, δ)
    D = length(xref)
    Aδ = Matrix(H.a(xref)) + δ * I
    Hxδ = (x, p) -> H.H_x(x, p)
    Hpδ = (x, p) -> H.H_p(x, p) .+ δ .* p
    return FreidlinWentzellHamiltonian{false, D}(
        Hxδ, Hpδ; a = Base.Returns(Aδ), x_ref = collect(xref),
    )
end

function _regularized_sgmam_seed(H, x_init)
    path = Matrix(x_init)
    result = nothing
    xref = collect(view(path, :, 1))
    optimizer = AdaptiveGeometricGradient(; stepsize = 100.0, probe_length = 25)
    for δ in (0.3, 0.1, 0.03)
        Hδ = _regularized_hamiltonian(H, xref, δ)
        result = minimize_geometric_action(
            Hδ, path, optimizer; maxiters = 200, show_progress = false,
        )
        path = CT._path_matrix(result.path)
    end
    return result
end

function _project_phase_endpoint(lin, x, p, eps_lin::T) where {T}
    y_raw = vcat(T.(x), T.(p)) .- lin.xstar_aug
    c_raw = lin.U' * y_raw
    norm_lin = LinearAlgebra.norm(lin.U * c_raw)
    scale = norm_lin > eps(T) ? eps_lin / norm_lin : eps_lin
    return T.(c_raw .* scale)
end

# Diagnostic override: use the global regularized-sgMAM continuation only to initialize the
# exact singular Hamiltonian BVP. No regularization enters the shooting residual or action.
@eval CT function _initial_guess_unknowns(ws::MultipleShootingWorkspace{IIP, D}, x_init) where {IIP, D}
    T = eltype(ws.grid)
    if !_rank_deficient_at_reference(ws, T.(x_init[:, 1]))
        N = size(x_init, 2)
        x_a_near = T.(x_init[:, min(2, N)])
        x_b_near = T.(x_init[:, max(N - 1, 1)])
        L0 = _initial_path_length(x_init, T)
        c_a = _project_endpoint(ws.lin_a, x_a_near, ws.eps_lin)
        c_b = _project_endpoint(ws.lin_b, x_b_near, ws.eps_lin)
        return _initial_guess_full_rank(ws, x_init, c_a, c_b, L0)
    end

    seed = Main._regularized_sgmam_seed(ws.H, x_init)
    path = _path_matrix(seed.path)
    p = seed.generalized_momentum
    N = size(path, 2)
    c_a = Main._project_phase_endpoint(ws.lin_a, view(path, :, min(2, N)), view(p, :, min(2, N)), ws.eps_lin)
    c_b = Main._project_phase_endpoint(ws.lin_b, view(path, :, max(N - 1, 1)), view(p, :, max(N - 1, 1)), ws.eps_lin)
    L0 = _initial_path_length(path, T)

    interior = zeros(T, 2D * (ws.nshoots - 1))
    for i in 1:(ws.nshoots - 1)
        idx = clamp(round(Int, (i / ws.nshoots) * (N - 1)) + 1, 1, N)
        @inbounds for k in 1:D
            interior[(i - 1) * 2D + k] = path[k, idx]
            interior[(i - 1) * 2D + D + k] = p[k, idx]
        end
    end
    return vcat(c_a, interior, c_b, [T(max(L0, ws.eps_lin))])
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

    @test_throws ArgumentError minimize_geometric_action(
        H, init, GeometricGradient(); maxiters = 1, show_progress = false
    )

    res = minimize_geometric_action(
        H,
        init,
        MultipleShooting(
            ; nshoots = 10, nlsolve = _rank_deficient_trust_region(), maxiters = 200,
            eps_lin = 1.0e-6, abstol = 1.0e-8, reltol = 1.0e-7,
        );
        show_progress = false,
    )

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
            ; nshoots = 10, nlsolve = _rank_deficient_trust_region(), maxiters = 200,
            eps_lin = 1.0e-6, abstol = 1.0e-8, reltol = 1.0e-7,
        );
        show_progress = false,
    )

    @test isapprox(res.action, 0.75; rtol = 3.0e-2)
    @test _H_invariant_max(H, res) < 1.0e-5
end
