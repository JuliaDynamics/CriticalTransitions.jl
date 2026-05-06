"""
Abstract supertype of boundary conditions accepted by [`DiffusionGenerator`](@ref).
Concrete subtypes: [`Reflecting`](@ref), [`Periodic`](@ref), [`Absorbing`](@ref).
"""
abstract type BoundaryCondition end

"""
No-flux Neumann boundary: outer faces are omitted from the stencil. The
chain bounces off the grid edges, mass is preserved.
"""
struct Reflecting <: BoundaryCondition end

"""
Periodic boundary: outer faces wrap to the opposite end of the same axis.
Use for systems on a torus / ring (angle variables).
"""
struct Periodic <: BoundaryCondition end

"""
Absorbing boundary: cells on the outer boundary leak probability through
their outer faces (escape rate from the cell-center drift). Row sums on
boundary cells become negative; the resulting `Q` is a sub-generator.
"""
struct Absorbing <: BoundaryCondition end

"""
$(TYPEDEF)

The discretised infinitesimal generator of an autonomous Itô diffusion on a
[`CartesianGrid`](@ref), built by the Scharfetter–Gummel exponential-
fitting finite-volume scheme.

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

Use [`m_matrix`](@ref) for the M-matrix view (PDE convention),
[`rate_matrix`](@ref) for the CTMC view, or
[`fokker_planck_operator`](@ref) for the forward Kolmogorov / FP operator
(`Qᵀ`).

# Fields
$(TYPEDFIELDS)

# Construction
```julia
DiffusionGenerator(sys::CoupledSDEs, grid::CartesianGrid; bc=Reflecting())
```
Discretises the SDE's backward Kolmogorov operator on `grid`. Drift `b(x)`
is read via [`drift`](@ref); diffusion is read once from
`covariance_matrix(sys)` and must be diagonal (axis-aligned noise).

Each axis-aligned face contributes

    Q[n, m] = ε / h_k² · B(-z) ,    Q[m, n] = ε / h_k² · B(z)

with `ε = Σ_kk / 2`, `z = b_k(x_face) · h_k / ε`, `B(z) = z / (eᶻ - 1)` the
Bernoulli function. The drift component `b_k` at the face center is averaged
from its two adjacent cell-center values.

`bc` controls how the outer grid faces are treated and is a singleton
[`BoundaryCondition`](@ref) instance applied to every axis, or a
`D`-tuple for per-axis control:

- [`Reflecting`](@ref)`()` (default) — outer faces are omitted; chain
  bounces off the edges.
- [`Periodic`](@ref)`()` — outer faces wrap to the opposite end of the
  same axis.
- [`Absorbing`](@ref)`()` — boundary cells leak probability through their
  outer faces. `Q` is a sub-generator (row sums on the boundary become
  negative). Useful for survival problems.
"""
struct DiffusionGenerator{D, BC <: Tuple}
    "The Cartesian grid the generator lives on."
    grid::CartesianGrid{D}
    "Sparse rate matrix (off-diagonals = transition rates ≥ 0)."
    Q::SparseMatrixCSC{Float64, Int}
    "Per-axis boundary condition (one [`BoundaryCondition`](@ref) per axis)."
    bc::BC
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
    tol = 1.0e-10 * max(maximum(abs, diag(Σ)), 1.0)
    @inbounds for i in 1:D, j in 1:D
        if i != j && abs(Σ[i, j]) > tol
            throw(
                ArgumentError(
                    "DiffusionGenerator requires diagonal noise covariance; got Σ[$i,$j] = $(Σ[i, j])",
                ),
            )
        end
    end
    return SVector{D, Float64}(ntuple(k -> Σ[k, k], D))
end

function _normalize_bc(bc::BoundaryCondition, ::Val{D}) where {D}
    return ntuple(_ -> bc, D)
end
function _normalize_bc(bc::Tuple, ::Val{D}) where {D}
    length(bc) == D || throw(
        ArgumentError(
            "expected a single BoundaryCondition or a $D-tuple thereof; " *
                "got a $(length(bc))-tuple",
        ),
    )
    @inbounds for k in 1:D
        bc[k] isa BoundaryCondition || throw(
            ArgumentError(
                "bc[$k] must be a BoundaryCondition (Reflecting/Periodic/" *
                    "Absorbing); got $(typeof(bc[k]))",
            ),
        )
    end
    return bc
end
_normalize_bc(bc, ::Val{D}) where {D} = throw(
    ArgumentError(
        "bc must be a BoundaryCondition or a $D-tuple of them; got $(typeof(bc))",
    ),
)

# Internal builder for the SG finite-volume sparse matrix.
# `sign = +1` produces the rate matrix / generator (off-diagonals ≥ 0).
# `sign = -1` produces the M-matrix form (off-diagonals ≤ 0).
function _assemble_generator(
        sys::CoupledSDEs,
        grid::CartesianGrid{D},
        sign::Int,
        bc::Tuple,
    )::SparseMatrixCSC{Float64, Int} where {D}
    sign == +1 || sign == -1 || throw(ArgumentError("sign must be ±1"))
    diffusion = _diagonal_diffusion(sys, Val(D))
    nbox = grid.nbox
    N = ncells(grid)
    LI = LinearIndices(nbox)

    drift_components = ntuple(_ -> Array{Float64, D}(undef, nbox), D)
    @inbounds for I in CartesianIndices(nbox)
        x = cell_center(grid, I)
        f = drift(sys, x)
        for k in 1:D
            drift_components[k][I] = f[k]
        end
    end

    # Generous upper bound on non-zeros: 2 directed edges per face + 1 diagonal
    # per cell. With periodic, +k boundary cells get one extra wraparound face
    # along that axis.
    nz_max = 2 * D * N + N
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
        bc_k = bc[k]

        # Helper closure for the SG (or upwind) face rates given a face drift.
        face_rates(bf::Float64) = if Dk > 0
            z = bf * h_ε
            (ε_h2 * _bernoulli(-z), ε_h2 * _bernoulli(z))
        else
            (max(bf, 0.0) / hk, max(-bf, 0.0) / hk)
        end

        for I in CartesianIndices(nbox)
            if I[k] < nbox[k]
                # Standard interior face.
                J = CartesianIndex(ntuple(d -> d == k ? I[d] + 1 : I[d], D))
                n = LI[I]
                m = LI[J]
                rate_nm, rate_mn = face_rates(0.5 * (fk[I] + fk[J]))

                idx += 1
                rows[idx] = n; cols[idx] = m; vals[idx] = sign * rate_nm
                idx += 1
                rows[idx] = m; cols[idx] = n; vals[idx] = sign * rate_mn
                diagacc[n] -= sign * rate_nm
                diagacc[m] -= sign * rate_mn
            else
                # I[k] == nbox[k]: handle boundary on the +k side.
                if bc_k isa Reflecting
                    # nothing to do
                elseif bc_k isa Periodic
                    # Wraparound face from the last cell along axis k to the
                    # first. Drift averaged across the wrap (treated like an
                    # ordinary face).
                    J = CartesianIndex(ntuple(d -> d == k ? 1 : I[d], D))
                    n = LI[I]
                    m = LI[J]
                    rate_nm, rate_mn = face_rates(0.5 * (fk[I] + fk[J]))

                    idx += 1
                    rows[idx] = n; cols[idx] = m; vals[idx] = sign * rate_nm
                    idx += 1
                    rows[idx] = m; cols[idx] = n; vals[idx] = sign * rate_mn
                    diagacc[n] -= sign * rate_nm
                    diagacc[m] -= sign * rate_mn
                elseif bc_k isa Absorbing
                    # Outflow through the +k wall using the cell-center drift
                    # as a face-drift approximation. Mass leaves the system.
                    n = LI[I]
                    rate_out, _ = face_rates(fk[I])
                    diagacc[n] -= sign * rate_out
                end
            end
        end

        # -k wall: only matters for absorbing (reflecting omits it; periodic
        # is already handled by the wraparound face above).
        if bc_k isa Absorbing
            for I in CartesianIndices(nbox)
                I[k] == 1 || continue
                n = LI[I]
                _, rate_out = face_rates(fk[I])
                diagacc[n] -= sign * rate_out
            end
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

function DiffusionGenerator(
        sys::CoupledSDEs, grid::CartesianGrid{D}; bc = Reflecting()
    ) where {D}
    bc_tuple = _normalize_bc(bc, Val(D))
    Q = _assemble_generator(sys, grid, +1, bc_tuple)
    return DiffusionGenerator{D, typeof(bc_tuple)}(grid, Q, bc_tuple)
end

"""
$(TYPEDSIGNATURES)

Rate matrix of `gen`: returns `gen.Q`. The same matrix is the **discrete
backward Kolmogorov / generator** — it acts on observables `u` via
`du/dt = Q u`, and the committor / mean first-passage time / other
backward-Kolmogorov BVPs are linear systems built on it.

In CTMC terminology the off-diagonals are transition rates between
adjacent cells, the diagonal is minus the escape rate, and rows sum to
zero. For the PDE M-matrix sign convention see [`m_matrix`](@ref); for the
forward Kolmogorov / Fokker–Planck operator (acts on densities), see
[`fokker_planck_operator`](@ref).
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

"""
$(TYPEDSIGNATURES)

Discrete **Fokker–Planck operator** (forward Kolmogorov operator) of
`gen`: returns `transpose(gen.Q)` materialised as a fresh sparse matrix.

It acts on probability densities: the discrete Fokker–Planck equation is

    dρ/dt = fokker_planck_operator(gen) * ρ ,

and the invariant density solves `fokker_planck_operator(gen) * ρ = 0`
(plus the normalisation `dot(ρ, fill(cell_volume(gen), N)) = 1`). For the
same density `ρ` and the relevant boundary conditions, the dual quantities
on observables (committors, MFPT) come from the backward operator
[`rate_matrix`](@ref).
"""
fokker_planck_operator(gen::DiffusionGenerator) = sparse(transpose(gen.Q))
