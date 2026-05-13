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

The discretised generator of an autonomous Itô diffusion on a
[`CartesianGrid`](@ref), built by the Scharfetter–Gummel exponential-
fitting finite-volume scheme.

The matrix `Q` generates the **transition semigroup** `exp(t Q)` of the
discretised process. With probability densities and observables represented
as column vectors (the convention used throughout this module), it
propagates observables forward (`du/dt = Q u`, the discrete backward
Kolmogorov) and densities forward via the transpose (`dρ/dt = Qᵀ ρ`, the
discrete Fokker–Planck, exposed as [`fokker_planck_operator`](@ref)).
Equivalently, `Q[i, j]` for `i ≠ j` is the transition **rate** from cell
`i` to cell `j`.

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


`bc` controls how the outer grid faces are treated and is a singleton
[`BoundaryCondition`](@ref) instance applied to every axis, or a
`D`-tuple for per-axis control:

- [`Reflecting`](@ref)`()` (default) — outer faces are omitted; chain
  bounces off the edges.
- [`Periodic`](@ref)`()` — outer faces wrap to the opposite end of the
  same axis.
- [`Absorbing`](@ref)`()` — boundary cells leak probability through their
  outer faces. `Q` is a sub-generator (row sums on the boundary become
  negative).
"""
struct DiffusionGenerator{D, BC <: Tuple, T <: AbstractFloat}
    "The Cartesian grid the generator lives on."
    grid::CartesianGrid{D, T}
    "Sparse rate matrix (off-diagonals = transition rates ≥ 0)."
    Q::SparseMatrixCSC{T, Int}
    "Per-axis boundary condition (one [`BoundaryCondition`](@ref) per axis)."
    bc::BC
end

@inline ncells(gen::DiffusionGenerator) = ncells(gen.grid)
@inline cell_volume(gen::DiffusionGenerator) = cell_volume(gen.grid)
@inline floattype(::DiffusionGenerator{D, BC, T}) where {D, BC, T} = T
@inline floattype(::CartesianGrid{D, T}) where {D, T} = T

function _diagonal_diffusion(sys::CoupledSDEs, ::Val{D}, ::Type{T}) where {D, T}
    Σ = covariance_matrix(sys)
    size(Σ) == (D, D) || throw(
        DimensionMismatch(
            "covariance matrix size $(size(Σ)) does not match grid dimension $D"
        ),
    )
    tol = T(1.0e-10) * max(maximum(abs, diag(Σ)), one(T))
    @inbounds for i in 1:D, j in 1:D
        if i != j && abs(Σ[i, j]) > tol
            throw(
                ArgumentError(
                    "DiffusionGenerator requires diagonal noise covariance; got Σ[$i,$j] = $(Σ[i, j])",
                ),
            )
        end
    end
    @inbounds for k in 1:D
        Σ[k, k] >= 0 || throw(
            ArgumentError("diffusion covariance diagonal Σ[$k,$k] must be non-negative"),
        )
    end
    return SVector{D, T}(ntuple(k -> T(Σ[k, k]), D))
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

# Assemble the SG finite-volume sparse matrix. `sign = -1` gives the
# corresponding M-matrix convention.
function _assemble_generator(
        sys::CoupledSDEs,
        grid::CartesianGrid{D, T},
        sign::Int,
        bc::BC,
    )::SparseMatrixCSC{T, Int} where {D, T, BC <: Tuple}
    sign == +1 || sign == -1 || throw(ArgumentError("sign must be ±1"))
    diffusion = _diagonal_diffusion(sys, Val(D), T)
    nbox = grid.nbox
    N = ncells(grid)
    LI = LinearIndices(nbox)

    drift_components = ntuple(_ -> Array{T, D}(undef, nbox), Val(D))
    @inbounds for I in CartesianIndices(nbox)
        x = cell_center(grid, I)
        f = drift(sys, x)
        for k in 1:D
            drift_components[k][I] = T(f[k])
        end
    end

    # Upper bound: two directed entries per face plus one diagonal per cell.
    nz_max = 2 * D * N + N
    rows = Vector{Int}(undef, nz_max)
    cols = Vector{Int}(undef, nz_max)
    vals = Vector{T}(undef, nz_max)
    diagacc = zeros(T, N)
    idx = 0

    # Function barrier: each call specialises on the concrete type of `bc[k]`,
    # turning the per-cell BC isa checks into compile-time branches.
    for k in 1:D
        idx = _fill_axis!(
            bc[k], rows, cols, vals, diagacc, idx, k, sign,
            drift_components[k], diffusion[k], grid.h[k], nbox, LI,
        )
    end

    @inbounds for n in 1:N
        idx += 1
        rows[idx] = n
        cols[idx] = n
        vals[idx] = diagacc[n]
    end

    return sparse(view(rows, 1:idx), view(cols, 1:idx), view(vals, 1:idx), N, N)
end

# SG rates for diffusive faces; upwind rates for deterministic axes.
@inline function _face_rates(bf::T, Dk::T, ε_h2::T, h_ε::T, hk::T) where {T}
    return if Dk > 0
        z = bf * h_ε
        (ε_h2 * _bernoulli(-z), ε_h2 * _bernoulli(z))
    else
        (max(bf, zero(T)) / hk, max(-bf, zero(T)) / hk)
    end
end

@inline function _fill_axis!(
        bc_k::BC_K,
        rows::Vector{Int}, cols::Vector{Int}, vals::Vector{T},
        diagacc::Vector{T}, idx::Int, k::Int, sign::Int,
        fk::AbstractArray{T, D}, Dk::T, hk::T,
        nbox::NTuple{D, Int}, LI::LinearIndices{D},
    ) where {T, D, BC_K <: BoundaryCondition}
    ε = Dk / 2
    ε_h2 = ε / (hk * hk)
    h_ε = ε > 0 ? hk / ε : T(NaN)

    @inbounds for I in CartesianIndices(nbox)
        if I[k] < nbox[k]
            J = CartesianIndex(Base.setindex(I.I, I[k] + 1, k))
            n = LI[I]
            m = LI[J]
            rate_nm, rate_mn = _face_rates(T(1) / 2 * (fk[I] + fk[J]), Dk, ε_h2, h_ε, hk)

            idx += 1
            rows[idx] = n; cols[idx] = m; vals[idx] = sign * rate_nm
            idx += 1
            rows[idx] = m; cols[idx] = n; vals[idx] = sign * rate_mn
            diagacc[n] -= sign * rate_nm
            diagacc[m] -= sign * rate_mn
        else
            if bc_k isa Periodic
                J = CartesianIndex(Base.setindex(I.I, 1, k))
                n = LI[I]
                m = LI[J]
                rate_nm, rate_mn = _face_rates(T(1) / 2 * (fk[I] + fk[J]), Dk, ε_h2, h_ε, hk)

                idx += 1
                rows[idx] = n; cols[idx] = m; vals[idx] = sign * rate_nm
                idx += 1
                rows[idx] = m; cols[idx] = n; vals[idx] = sign * rate_mn
                diagacc[n] -= sign * rate_nm
                diagacc[m] -= sign * rate_mn
            elseif bc_k isa Absorbing
                n = LI[I]
                rate_out, _ = _face_rates(fk[I], Dk, ε_h2, h_ε, hk)
                diagacc[n] -= sign * rate_out
            end
            # Reflecting: nothing to do.
        end
    end

    # -k wall only contributes for absorbing boundaries.
    if bc_k isa Absorbing
        @inbounds for I in CartesianIndices(nbox)
            I[k] == 1 || continue
            n = LI[I]
            _, rate_out = _face_rates(fk[I], Dk, ε_h2, h_ε, hk)
            diagacc[n] -= sign * rate_out
        end
    end

    return idx
end

function DiffusionGenerator(
        sys::CoupledSDEs, grid::CartesianGrid{D, T}; bc = Reflecting()
    ) where {D, T}
    bc_tuple = _normalize_bc(bc, Val(D))
    Q = _assemble_generator(sys, grid, +1, bc_tuple)
    return DiffusionGenerator{D, typeof(bc_tuple), T}(grid, Q, bc_tuple)
end

"""
$(TYPEDSIGNATURES)

Rate matrix of `gen`: returns `gen.Q`, also called the **discrete
backward Kolmogorov / generator**.

The off-diagonals are transition rates between adjacent cells, the diagonal is minus the
escape rate. For the PDE M-matrix sign convention see [`m_matrix`](@ref); for the
forward Kolmogorov / Fokker–Planck operator (acts on densities), see
[`fokker_planck_operator`](@ref).
"""
@inline rate_matrix(generator::DiffusionGenerator) = generator.Q

"""
$(TYPEDSIGNATURES)

M-matrix of the `generator`. Equivalent to `-rate_matrix(gen)`.
"""
m_matrix(generator::DiffusionGenerator) = -generator.Q

"""
$(TYPEDSIGNATURES)

Discrete **Fokker–Planck operator** (forward Kolmogorov operator) of
`generator`.

It acts on probability densities: the discrete Fokker–Planck equation is

    dρ/dt = fokker_planck_operator(gen) * ρ ,

and the invariant density solves `fokker_planck_operator(gen) * ρ = 0`.
"""
fokker_planck_operator(generator::DiffusionGenerator) = sparse(transpose(generator.Q))
