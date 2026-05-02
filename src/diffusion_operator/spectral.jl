# =====================================================================
# Spectral analysis of a DiffusionGenerator. The eigenvalues of `Q` have
# non-positive real part; the largest is exactly zero with eigenvector
# `1`, and the next ones decay at rates set by the spectral gap.
# =====================================================================

"""
$(TYPEDSIGNATURES)

Slowest-decaying eigenmodes of the generator `gen`. Returns
`(λ, V)` where `λ::Vector{ComplexF64}` are the `k` eigenvalues with
largest real part (least negative — slowest decay), sorted descending,
and `V::Matrix{ComplexF64}` collects the corresponding right eigenvectors
of `Q` columnwise (`Q * V[:, i] = λ[i] * V[:, i]`).

The eigenvalue closest to zero is exactly `0`, with eigenvector the
constant function (since `Q * 1 = 0`); the corresponding *left*
eigenvector is the invariant density. The remaining eigenvalues lie in
the open left half-plane, and `-1/Re(λ_k)` gives the metastable
timescale of the `k`-th mode. The slow non-trivial right eigenvectors
are the canonical reaction coordinates / metastable-state indicators
used in Markov state model decomposition (e.g. PCCA+).

Currently uses dense [`LinearAlgebra.eigen`](@extref) on `Matrix(gen.Q)`.
This is fine up to `ncells(gen)` of a few thousand; for larger grids,
run an iterative sparse eigensolver (Arpack.jl / KrylovKit.jl) directly
on [`rate_matrix(gen)`](@ref) — the slowest modes correspond to
`which = :LR` (largest real part).
"""
function eigenmodes(gen::DiffusionGenerator, k::Integer = 10)
    N = size(gen.Q, 1)
    k >= 1 || throw(ArgumentError("k must be ≥ 1"))
    k = min(k, N)
    F = eigen(Matrix(gen.Q))
    perm = sortperm(F.values; by = real, rev = true)
    inds = perm[1:k]
    return F.values[inds], F.vectors[:, inds]
end
