"""
$(TYPEDEF)

The discretised infinitesimal generator of an autonomous Itô diffusion on a
[`CartesianGrid`](@ref). Constructed by the Scharfetter–Gummel exponential-
fitting finite-volume scheme with reflected (no-flux Neumann) boundary
conditions.

The matrix `Q` generates the **transition semigroup** `exp(t Q)` of the
discretised process: it propagates probability distributions in time
(`dρ/dt = ρ Q`, the discrete Fokker–Planck) and observables (`du/dt = Q u`,
the discrete backward Kolmogorov). Equivalently, `Q[i, j]` for `i ≠ j` is
the transition **rate** from cell `i` to cell `j`.

The same matrix has several traditional names:

| Name | Where it comes from | Sign convention of `Q` |
|---|---|---|
| **rate matrix** / **Q-matrix** | probability theory, Markov chains | off-diagonals ≥ 0 |
| **infinitesimal generator** | semigroup theory, applied math | off-diagonals ≥ 0 |
| **negative M-matrix** | numerical PDE / finite-volume / FEM | (the M-matrix is `-Q`, off-diag ≤ 0) |

Use [`m_matrix`](@ref) to obtain the M-matrix view (positive diagonal,
non-positive off-diagonals) when the PDE/FV convention is more natural,
or use [`rate_matrix`](@ref) (alias for the field) when emphasising the
CTMC interpretation. They refer to the same operator with opposite signs.

# Fields
$(TYPEDFIELDS)

# Construction
```julia
DiffusionGenerator(sys::CoupledSDEs, grid::CartesianGrid)
```
Discretises the SDE's backward Kolmogorov operator on `grid`. Drift `b(x)`
is read via [`drift`](@ref); diffusion is read once from
`covariance_matrix(sys)` and must be diagonal (axis-aligned noise).

Each axis-aligned face contributes

    Q[n, m] = ε / h_k² · B(-z) ,    Q[m, n] = ε / h_k² · B(z)

with `ε = Σ_kk / 2`, `z = b_k(x_face) · h_k / ε`, `B(z) = z / (eᶻ - 1)` the
Bernoulli function. The drift component `b_k` at the face center is averaged
from its two adjacent cell-center values. Faces at the outer grid boundary
are simply omitted.
"""
struct DiffusionGenerator{D}
    "The Cartesian grid the generator lives on."
    grid::CartesianGrid{D}
    "Sparse rate matrix (off-diagonals = transition rates ≥ 0, row sums zero)."
    Q::SparseMatrixCSC{Float64,Int}
end

@inline ncells(gen::DiffusionGenerator) = ncells(gen.grid)
@inline cell_volume(gen::DiffusionGenerator) = cell_volume(gen.grid)

function _diagonal_diffusion(sys::CoupledSDEs, ::Val{D}) where {D}
    Σ = covariance_matrix(sys)
    size(Σ) == (D, D) || throw(
        DimensionMismatch(
            "covariance matrix size $(size(Σ)) does not match grid dimension $D"
        ),
    )
    tol = 1e-10 * max(maximum(abs, diag(Σ)), 1.0)
    @inbounds for i in 1:D, j in 1:D
        if i != j && abs(Σ[i, j]) > tol
            throw(
                ArgumentError(
                    "DiffusionGenerator requires diagonal noise covariance; got Σ[$i,$j] = $(Σ[i,j])",
                ),
            )
        end
    end
    return SVector{D,Float64}(ntuple(k -> Σ[k, k], D))
end

# Internal builder for the SG finite-volume sparse matrix.
# `sign = +1` produces the rate matrix / generator (off-diagonals ≥ 0).
# `sign = -1` produces the M-matrix form (off-diagonals ≤ 0).
function _assemble_generator(
    sys::CoupledSDEs, grid::CartesianGrid{D}, sign::Int
)::SparseMatrixCSC{Float64,Int} where {D}
    sign == +1 || sign == -1 || throw(ArgumentError("sign must be ±1"))
    diffusion = _diagonal_diffusion(sys, Val(D))
    nbox = grid.nbox
    N = ncells(grid)
    LI = LinearIndices(nbox)

    drift_components = ntuple(_ -> Array{Float64,D}(undef, nbox), D)
    @inbounds for I in CartesianIndices(nbox)
        x = cell_center(grid, I)
        f = drift(sys, x)
        for k in 1:D
            drift_components[k][I] = f[k]
        end
    end

    nfaces = sum(N - N ÷ nbox[k] for k in 1:D)
    nz_max = 2 * nfaces + N
    rows = Vector{Int}(undef, nz_max)
    cols = Vector{Int}(undef, nz_max)
    vals = Vector{Float64}(undef, nz_max)
    diagacc = zeros(Float64, N)
    idx = 0

    @inbounds for k in 1:D
        Dk = diffusion[k]
        hk = grid.h[k]
        ε = Dk / 2
        ε_h2 = ε / (hk * hk)
        h_ε = ε > 0 ? hk / ε : NaN
        fk = drift_components[k]

        for I in CartesianIndices(nbox)
            I[k] == nbox[k] && continue
            J = CartesianIndex(ntuple(d -> d == k ? I[d] + 1 : I[d], D))
            n = LI[I]
            m = LI[J]

            b_face = 0.5 * (fk[I] + fk[J])

            if Dk > 0
                z = b_face * h_ε
                rate_nm = ε_h2 * _bernoulli(-z)
                rate_mn = ε_h2 * _bernoulli(z)
            else
                rate_nm = max(b_face, 0.0) / hk
                rate_mn = max(-b_face, 0.0) / hk
            end

            idx += 1
            rows[idx] = n
            cols[idx] = m
            vals[idx] = sign * rate_nm
            idx += 1
            rows[idx] = m
            cols[idx] = n
            vals[idx] = sign * rate_mn

            diagacc[n] -= sign * rate_nm
            diagacc[m] -= sign * rate_mn
        end
    end

    @inbounds for n in 1:N
        idx += 1
        rows[idx] = n
        cols[idx] = n
        vals[idx] = diagacc[n]
    end

    return sparse(view(rows, 1:idx), view(cols, 1:idx), view(vals, 1:idx), N, N)
end

function DiffusionGenerator(sys::CoupledSDEs, grid::CartesianGrid{D}) where {D}
    Q = _assemble_generator(sys, grid, +1)
    return DiffusionGenerator{D}(grid, Q)
end

"""
$(TYPEDSIGNATURES)

Rate matrix (CTMC convention) of `gen`: returns `gen.Q`. Alias for direct
field access, named to match probability-theory terminology — off-diagonals
are transition rates between adjacent cells, the diagonal is the negative
escape rate, rows sum to zero.

For the same operator in PDE / M-matrix sign convention, see
[`m_matrix`](@ref).
"""
@inline rate_matrix(gen::DiffusionGenerator) = gen.Q

"""
$(TYPEDSIGNATURES)

M-matrix view of `gen`: returns `-gen.Q`, a sparse matrix with positive
diagonal, non-positive off-diagonals, and row sums zero. The committor
satisfies the homogeneous linear system `(M-matrix) q = 0` with Dirichlet
boundary conditions; the same operator is used for mean first-passage
times, spectral analysis, and other backward-Kolmogorov problems.

Equivalent to `-rate_matrix(gen)`. Allocates a fresh sparse matrix.
"""
m_matrix(gen::DiffusionGenerator) = -gen.Q
