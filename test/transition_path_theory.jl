using Test
using SparseArrays: rowvals, nonzeros, nzrange

using CriticalTransitions

# =====================================================================
# Helpers
# =====================================================================

# Cross-check the cached rate against an independent sum over the B boundary:
# k_AB = ∑_{i ∉ B, j ∈ B} ρ[i] v[i] q⁻[i] Q[i, j] q⁺[j].
function _rate_into_B(res::ReactiveTransition)
    grid = res.generator.grid
    v = prod(grid.h)
    Q = res.generator.Q
    rv = rowvals(Q)
    nz = nonzeros(Q)
    k = 0.0
    for col in 1:size(Q, 2)
        res.B_mask[col] || continue
        for p in nzrange(Q, col)
            row = rv[p]
            row != col || continue
            res.B_mask[row] && continue
            k += res.ρ[row] * v * res.qminus[row] * nz[p] * res.qplus[col]
        end
    end
    return k
end

# =====================================================================
# Reactive rate self-consistency
# =====================================================================

@testset "Self-consistent reactive rate" begin
    sys = CoupledSDEs((u, p, t) -> [u[1] - u[1]^3], [0.0]; noise_strength=0.6)
    grid = CartesianGrid((-2.0, 2.0, 200))
    res = ReactiveTransition(sys, grid, x -> x[1] < -0.7, x -> x[1] > 0.7)
    k_outA = reactive_rate(res)
    k_intoB = _rate_into_B(res)
    @test k_outA > 0
    @test isfinite(k_outA)
    @test isapprox(k_outA, k_intoB; rtol=1e-6)
end

# =====================================================================
# Higher-dimensional dispatch
# =====================================================================

@testset "3D dispatch" begin
    sys = CoupledSDEs(
        (u, p, t) -> [-u[1], -u[2], -u[3]], [0.0, 0.0, 0.0]; noise_strength=1.0
    )
    grid = CartesianGrid((-3.0, 3.0, 21), (-3.0, 3.0, 21), (-3.0, 3.0, 21))
    res = ReactiveTransition(
        sys,
        grid,
        x -> sum((x .- (-1.0, 0.0, 0.0)) .^ 2) < 0.4^2,
        x -> sum((x .- (1.0, 0.0, 0.0)) .^ 2) < 0.4^2,
    )
    @test length(stationary_distribution(res)) == 21^3
    @test size(res.generator.Q) == (21^3, 21^3)
    @test maximum(abs, vec(sum(res.generator.Q; dims=2))) < 1e-10
    @test reactive_rate(res) > 0
    @test extrema(forward_committor(res)) == (0.0, 1.0)
    J_nodes, J_faces = reactive_current(res)
    @test length(J_nodes) == 3 && length(J_faces) == 3
    @test size(J_nodes[1]) == (21, 21, 21)
    @test size(J_faces[1]) == (20, 21, 21)
    @test size(J_faces[2]) == (21, 20, 21)
    @test size(J_faces[3]) == (21, 21, 20)
end

# =====================================================================
# Probability accessors
# =====================================================================

@testset "Probability accessors" begin
    sys = CoupledSDEs((u, p, t) -> [u[1] - u[1]^3], [0.0]; noise_strength=0.6)
    grid = CartesianGrid((-2.0, 2.0, 200))
    res = ReactiveTransition(sys, grid, x -> x[1] < -0.7, x -> x[1] > 0.7)
    @test probability_reactive(res) >= 0
    @test probability_last_A(res) ≈ 0.5 atol = 1e-2
end

# =====================================================================
# Generator reuse across (A, B) sweeps
# =====================================================================

@testset "Generator reuse" begin
    sys = CoupledSDEs((u, p, t) -> [-u[1]], [0.0]; noise_strength=1.0)
    grid = CartesianGrid((-4.0, 4.0, 200))
    gen = DiffusionGenerator(sys, grid)
    r1 = ReactiveTransition(gen, x -> x[1] < -2, x -> x[1] > 2)
    r2 = ReactiveTransition(gen, x -> x[1] < -1, x -> x[1] > 1)
    @test r1.generator === gen
    @test r2.generator === gen
    @test reactive_rate(r1) < reactive_rate(r2)
end

# =====================================================================
# Physical reverse drift
# =====================================================================

@testset "Physical reverse drift" begin
    sys = CoupledSDEs((u, p, t) -> [-u[1]], [0.0]; noise_strength=1.0)
    grid = CartesianGrid((-4.0, 4.0, 200))
    res = ReactiveTransition(sys, grid, x -> x[1] < -2, x -> x[1] > 2; reverse=sys)
    @test res.physical_reverse
    res_chain = ReactiveTransition(sys, grid, x -> x[1] < -2, x -> x[1] > 2)
    @test maximum(abs.(res.qminus .- res_chain.qminus)) < 1e-6
end

# =====================================================================
# Reactive current decomposition
# =====================================================================

@testset "reactive_current with Periodic BCs includes wraparound face" begin
    sys = CoupledSDEs((u, p, t) -> [0.0], [0.0]; noise_strength=0.5)
    grid = CartesianGrid((-pi, pi, 60))

    # A near -π and B near +π so the geometry forces flux through the
    # wraparound face on the ring.
    A = x -> -3.0 < x[1] < -2.7
    B = x -> 2.7 < x[1] < 3.0
    res = ReactiveTransition(sys, grid, A, B; bc=Periodic())

    J_nodes, J_faces = reactive_current(res)
    @test size(J_faces[1]) == (60,)        # nbox[1] (not nbox[1] - 1)
    @test size(J_nodes[1]) == (60,)
    # Wraparound face carries a non-trivial flux.
    @test abs(J_faces[1][end]) > 1e-3

    # Reflecting BC keeps the old shape (nbox[k] - 1 along axis k).
    res_ref = ReactiveTransition(sys, grid, A, B)
    _, J_faces_ref = reactive_current(res_ref)
    @test size(J_faces_ref[1]) == (59,)
end

@testset "ReactiveTransition rejects Absorbing BC" begin
    sys = CoupledSDEs((u, p, t) -> [-u[1]], [0.0]; noise_strength=1.0)
    grid = CartesianGrid((-3.0, 3.0, 30))
    gen_abs = DiffusionGenerator(sys, grid; bc=Absorbing())
    @test_throws ArgumentError ReactiveTransition(gen_abs, x -> x[1] < -1, x -> x[1] > 1)
    @test_throws ArgumentError ReactiveTransition(
        sys, grid, x -> x[1] < -1, x -> x[1] > 1; bc=Absorbing()
    )
end

@testset "reactive_current rev/irr decomposition" begin
    # Reversible system: irreversible part should vanish.
    sys_rev = CoupledSDEs(
        (u, p, t) -> [u[1] - u[1]^3, -u[2]], [0.0, 0.0]; noise_strength=0.4
    )
    grid = CartesianGrid((-1.6, 1.6, 41), (-1.0, 1.0, 27))
    A = x -> sum((x .- (-1.0, 0.0)) .^ 2) < 0.25^2
    B = x -> sum((x .- (1.0, 0.0)) .^ 2) < 0.25^2
    res = ReactiveTransition(sys_rev, grid, A, B)

    Jn_full, Jf_full = reactive_current(res)
    Jn_rev, Jf_rev = reactive_current_reversible(res)
    Jn_irr, Jf_irr = reactive_current_irreversible(res)

    # Additivity
    for k in 1:2
        @test maximum(abs.(Jf_full[k] .- Jf_rev[k] .- Jf_irr[k])) < 1e-12
        @test maximum(abs.(Jn_full[k] .- Jn_rev[k] .- Jn_irr[k])) < 1e-12
    end
    # For reversible system, irrev part is essentially zero
    scale = maximum(maximum(abs, Jf_full[k]) for k in 1:2)
    @test maximum(maximum(abs, Jf_irr[k]) for k in 1:2) < 1e-10 * scale

    # Non-equilibrium system: irrev part should be substantial.
    sys_nq = CoupledSDEs(
        (u, p, t) -> [u[1] - u[1]^3 - 10 * u[1] * u[2]^2, -(1 + u[1]^2) * u[2]],
        [0.0, 0.0]; noise_strength=0.4,
    )
    res_nq = ReactiveTransition(sys_nq, grid, A, B)
    _, Jf_full_nq = reactive_current(res_nq)
    _, Jf_irr_nq = reactive_current_irreversible(res_nq)
    scale_nq = maximum(maximum(abs, Jf_full_nq[k]) for k in 1:2)
    @test maximum(maximum(abs, Jf_irr_nq[k]) for k in 1:2) > 0.05 * scale_nq
end

# =====================================================================
# Krylov solver via LinearSolve
# =====================================================================

@testset "ReactiveTransition with Krylov solver" begin
    using LinearSolve: KrylovJL_GMRES
    sys = CoupledSDEs((u, p, t) -> [u[1] - u[1]^3], [0.0]; noise_strength=0.6)
    grid = CartesianGrid((-2.0, 2.0, 200))
    gen = DiffusionGenerator(sys, grid)
    A = x -> x[1] < -0.7
    B = x -> x[1] > 0.7

    res_lu = ReactiveTransition(gen, A, B)
    res_kr = ReactiveTransition(gen, A, B; alg=KrylovJL_GMRES())
    @test reactive_rate(res_kr) ≈ reactive_rate(res_lu) rtol = 1e-8
end

# =====================================================================
# Type stability
# =====================================================================

@testset "Type stability (ReactiveTransition)" begin
    sys = CoupledSDEs((u, p, t) -> [-u[1]], [0.0]; noise_strength=1.0)
    grid = CartesianGrid((-4.0, 4.0, 50))
    gen = DiffusionGenerator(sys, grid)
    res = ReactiveTransition(gen, x -> x[1] < -2, x -> x[1] > 2)
    @test (@inferred forward_committor(res)) isa Vector{Float64}
    @test (@inferred backward_committor(res)) isa Vector{Float64}
    @test (@inferred stationary_distribution(res)) isa Vector{Float64}
    @test (@inferred reactive_density(res)) isa Vector{Float64}
    @test (@inferred reactive_rate(res)) isa Float64
    @test (@inferred probability_reactive(res)) isa Float64
    @test (@inferred probability_last_A(res)) isa Float64
end
