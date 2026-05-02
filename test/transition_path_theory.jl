using Test
using LinearAlgebra: diag
using SparseArrays: rowvals, nonzeros, nzrange, SparseMatrixCSC

using CriticalTransitions

# =====================================================================
# Helpers
# =====================================================================

function _offdiag_extrema(A)
    rv = rowvals(A)
    nz = nonzeros(A)
    omin, omax = Inf, -Inf
    for col in 1:size(A, 2), p in nzrange(A, col)
        row = rv[p]
        row != col || continue
        omin = min(omin, nz[p])
        omax = max(omax, nz[p])
    end
    return omin, omax
end

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
# DiffusionGenerator: shape and sign conventions
# =====================================================================

@testset "DiffusionGenerator (1D OU): rate-matrix structure" begin
    sys = CoupledSDEs((u, p, t) -> [-u[1]], [0.0]; noise_strength=1.0)
    grid = CartesianGrid((-4.0, 4.0, 200))
    gen = DiffusionGenerator(sys, grid)
    @test gen isa DiffusionGenerator
    Q = gen.Q
    @test Q isa SparseMatrixCSC{Float64,Int}
    @test size(Q) == (200, 200)
    @test maximum(abs, vec(sum(Q; dims=2))) < 1e-12
    omin, omax = _offdiag_extrema(Q)
    @test omin >= 0           # rate matrix: off-diagonals = transition rates ≥ 0
    @test isfinite(omax)
    @test all(diag(Q) .<= 0)  # rate matrix: diagonal = -escape rate ≤ 0
end

@testset "rate_matrix and m_matrix accessors" begin
    sys = CoupledSDEs((u, p, t) -> [-u[1]], [0.0]; noise_strength=1.0)
    grid = CartesianGrid((-4.0, 4.0, 200))
    gen = DiffusionGenerator(sys, grid)
    Q = rate_matrix(gen)
    L = m_matrix(gen)
    @test Q === gen.Q                              # alias for the field
    @test maximum(abs, L + Q) < 1e-12              # M-matrix is the negation
    omin, omax = _offdiag_extrema(L)
    @test omax <= 0           # M-matrix: off-diagonals ≤ 0
    @test all(diag(L) .>= 0)  # M-matrix: diagonal ≥ 0
end

# =====================================================================
# Stationary density
# =====================================================================

@testset "Invariant density (1D OU)" begin
    sys = CoupledSDEs((u, p, t) -> [-u[1]], [0.0]; noise_strength=1.0)
    grid = CartesianGrid((-4.0, 4.0, 200))
    res = ReactiveTransition(sys, grid, x -> x[1] < -2.0, x -> x[1] > 2.0)
    @test res isa ReactiveTransition
    ρ = stationary_distribution(res)
    @test sum(ρ) * grid.h[1] ≈ 1.0 atol = 1e-12
    xs = collect(grid.centers[1])
    ρ_analytic = exp.(-xs .^ 2) ./ sqrt(pi)
    @test sum(abs.(ρ .- ρ_analytic)) * grid.h[1] < 5e-3
end

@testset "Invariant density (2D OU)" begin
    sys = CoupledSDEs((u, p, t) -> [-u[1], -u[2]], [0.0, 0.0]; noise_strength=1.0)
    grid = CartesianGrid((-4.0, 4.0, 81), (-4.0, 4.0, 81))
    res = ReactiveTransition(sys, grid, x -> x[1] < -2.0, x -> x[1] > 2.0)
    ρ = stationary_distribution(res)
    @test sum(ρ) * prod(grid.h) ≈ 1.0 atol = 1e-12
    xs = [grid.centers[1][I[1]] for I in CartesianIndices(grid.nbox)]
    ys = [grid.centers[2][I[2]] for I in CartesianIndices(grid.nbox)]
    ρ_analytic = vec(exp.(-(xs .^ 2 .+ ys .^ 2)) ./ pi)
    @test sum(abs.(ρ .- ρ_analytic)) * prod(grid.h) < 5e-3
end

# =====================================================================
# Committor: BCs and reversibility
# =====================================================================

@testset "Committor BCs and reversibility (1D double-well)" begin
    sys = CoupledSDEs((u, p, t) -> [u[1] - u[1]^3], [0.0]; noise_strength=0.7)
    grid = CartesianGrid((-2.0, 2.0, 200))
    res = ReactiveTransition(sys, grid, x -> x[1] < -0.7, x -> x[1] > 0.7)
    qp = forward_committor(res)
    qm = backward_committor(res)
    @test all(qp[res.A_mask] .== 0)
    @test all(qp[res.B_mask] .== 1)
    @test all(qm[res.A_mask] .== 1)
    @test all(qm[res.B_mask] .== 0)
    @test extrema(qp) == (0.0, 1.0)
    @test extrema(qm) == (0.0, 1.0)
    @test maximum(abs.(qm .- (1 .- qp))) < 1e-6
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
# Layer-3 analyses on a DiffusionGenerator
# =====================================================================

@testset "Layer-3 analyses" begin
    sys = CoupledSDEs((u, p, t) -> [u[1] - u[1]^3], [0.0]; noise_strength=0.6)
    grid = CartesianGrid((-2.0, 2.0, 200))
    gen = DiffusionGenerator(sys, grid)

    @testset "stationary_distribution(gen) matches ReactiveTransition.ρ" begin
        ρ_gen = stationary_distribution(gen)
        res = ReactiveTransition(gen, x -> x[1] < -0.7, x -> x[1] > 0.7)
        @test maximum(abs.(ρ_gen .- res.ρ)) < 1e-12
    end

    @testset "forward_committor(gen, A, B) matches ReactiveTransition.qplus" begin
        A = x -> x[1] < -0.7
        B = x -> x[1] > 0.7
        qp = forward_committor(gen, A, B)
        res = ReactiveTransition(gen, A, B)
        @test maximum(abs.(qp .- res.qplus)) < 1e-12
        @test extrema(qp) == (0.0, 1.0)
    end

    @testset "backward_committor(gen, A, B) matches ReactiveTransition.qminus" begin
        A = x -> x[1] < -0.7
        B = x -> x[1] > 0.7
        qm = backward_committor(gen, A, B)
        res = ReactiveTransition(gen, A, B)
        @test maximum(abs.(qm .- res.qminus)) < 1e-12
    end

    @testset "backward_committor with explicit reverse generator" begin
        A = x -> x[1] < -0.7
        B = x -> x[1] > 0.7
        qm_explicit = backward_committor(gen, A, B; reverse=gen)
        qm_adjoint = backward_committor(gen, A, B)
        @test maximum(abs.(qm_explicit .- qm_adjoint)) < 1e-6
    end

    @testset "mean_first_passage_time" begin
        sys_ou = CoupledSDEs((u, p, t) -> [-u[1]], [0.0]; noise_strength=1.0)
        grid_ou = CartesianGrid((-4.0, 4.0, 200))
        gen_ou = DiffusionGenerator(sys_ou, grid_ou)
        τ = mean_first_passage_time(gen_ou, x -> abs(x[1]) > 2)
        @test all(τ .>= 0)
        target_mask = [abs(x) > 2 for x in grid_ou.centers[1]]
        @test maximum(abs.(τ[target_mask])) < 1e-10
        @test minimum(τ[.!target_mask]) > 0
        i0 = findfirst(x -> x >= 0, grid_ou.centers[1])
        @test τ[i0] ≈ maximum(τ) atol = 1e-3 * maximum(τ)
    end
end

# =====================================================================
# Type stability
# =====================================================================

@testset "Type stability" begin
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
    @test (@inferred stationary_distribution(gen)) isa Vector{Float64}
    @test (@inferred mean_first_passage_time(gen, x -> abs(x[1]) > 2)) isa Vector{Float64}
end
