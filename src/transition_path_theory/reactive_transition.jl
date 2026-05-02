"""
$(TYPEDEF)

Bundled transition-path-theory result on a [`DiffusionGenerator`](@ref).
Carries the generator, the metastable-set masks, and the cached committor /
invariant-density / rate solutions. The reactive rate is computed once at
construction and cached as a field; vector-valued derived quantities
(`reactive_density`, `reactive_current`) stay on demand.

Internally calls the Layer-3 analyses
([`stationary_distribution`](@ref), [`forward_committor`](@ref),
[`backward_committor`](@ref)) and bundles their outputs together with the
A → B reactive rate.

# Fields
$(TYPEDFIELDS)

# Construction
```julia
ReactiveTransition(gen::DiffusionGenerator, A, B; reverse=nothing)
ReactiveTransition(sys::CoupledSDEs, grid::CartesianGrid, A, B; reverse=nothing)
```
The first form takes a pre-built [`DiffusionGenerator`](@ref) — useful for
sweeping multiple `(A, B)` pairs without rediscretising. The second is a
convenience that builds the generator from `sys` internally.

`A` and `B` accept:
- a predicate `x::SVector{D} -> Bool` evaluated at each cell center,
- a `BitVector` of length `ncells(grid)`, or
- a vector of linear cell indices.

If `reverse` is provided (a `CoupledSDEs` for the high-level form, a
`DiffusionGenerator` for the generator form), the backward committor is
computed on that physical reverse-drift generator instead of the discrete
adjoint.

`alg = nothing` (default) uses sparse direct LU for all internal solves.
Pass a LinearSolve.jl algorithm (e.g. `KrylovJL_GMRES()`) to use an
iterative solver instead — useful for large grids where direct LU runs
out of memory.
"""
struct ReactiveTransition{D}
    "The diffusion generator the analysis was built on."
    generator::DiffusionGenerator{D}
    "Boolean mask selecting cells in set A (length `ncells(grid)`)."
    A_mask::BitVector
    "Boolean mask selecting cells in set B (length `ncells(grid)`)."
    B_mask::BitVector
    "Invariant probability density (`dot(ρ, weights) = 1`)."
    ρ::Vector{Float64}
    "Forward committor `q⁺` (probability of reaching B before A)."
    qplus::Vector{Float64}
    "Backward committor `q⁻`."
    qminus::Vector{Float64}
    "A → B reactive transition rate (cached)."
    rate::Float64
    "True if `qminus` came from a user-supplied physical reverse generator; false if from the discrete adjoint."
    physical_reverse::Bool
end

function ReactiveTransition(
        gen::DiffusionGenerator{D},
        A,
        B;
        reverse::Union{Nothing, DiffusionGenerator{D}} = nothing,
        alg = nothing,
    ) where {D}
    any(b -> b isa Absorbing, gen.bc) && throw(
        ArgumentError(
            "ReactiveTransition is undefined for absorbing boundary " *
                "conditions: the chain has mass loss and admits no nontrivial " *
                "invariant density. Use Reflecting() or Periodic() instead.",
        ),
    )
    grid = gen.grid
    A_mask = _to_mask(A, grid)
    B_mask = _to_mask(B, grid)
    any(A_mask) || throw(ArgumentError("set A is empty"))
    any(B_mask) || throw(ArgumentError("set B is empty"))
    any(A_mask .& B_mask) && throw(ArgumentError("sets A and B overlap"))

    Q = gen.Q
    weights = fill(cell_volume(grid), ncells(gen))

    ρ = _invariant_density(Q, weights; alg = alg)
    qplus = _forward_committor(Q, A_mask, B_mask; alg = alg)

    if reverse === nothing
        qminus = _backward_committor_adjoint(Q, ρ, A_mask, B_mask; alg = alg)
        physical_reverse = false
    else
        qminus = _backward_committor_explicit(reverse.Q, A_mask, B_mask; alg = alg)
        physical_reverse = true
    end

    rate = _reactive_rate(Q, ρ, qplus, qminus, weights, A_mask)
    return ReactiveTransition{D}(gen, A_mask, B_mask, ρ, qplus, qminus, rate, physical_reverse)
end

function ReactiveTransition(
        sys::CoupledSDEs,
        grid::CartesianGrid{D},
        A,
        B;
        reverse::Union{Nothing, CoupledSDEs} = nothing,
        bc = Reflecting(),
        alg = nothing,
    ) where {D}
    gen = DiffusionGenerator(sys, grid; bc = bc)
    gen_rev =
        reverse === nothing ? nothing : DiffusionGenerator(reverse, grid; bc = bc)
    return ReactiveTransition(gen, A, B; reverse = gen_rev, alg = alg)
end

# =====================================================================
# Cached accessors (constant-time field reads)
# =====================================================================

@inline forward_committor(r::ReactiveTransition) = r.qplus
@inline backward_committor(r::ReactiveTransition) = r.qminus
@inline stationary_distribution(r::ReactiveTransition) = r.ρ

"""
$(TYPEDSIGNATURES)

A → B **reactive transition rate** ``k_{AB}``: events per unit time, the
rate at which probability leaves `A` along reactive trajectories. Cached
on the [`ReactiveTransition`](@ref) — constant-time field read.

In CTMC form,

    k_{AB} = ∑_{i ∈ A, j ∉ A} ρ[i] · v · q⁻[i] · Q[i, j] · q⁺[j] ,

with `v = cell_volume(grid)`. By Vanden-Eijnden's identity the sum over
`{i ∉ B, j ∈ B}` (the boundary of `B`) gives the same value to solver
tolerance.
"""
@inline reactive_rate(r::ReactiveTransition) = r.rate

"""
$(TYPEDSIGNATURES)

Per-cell reactive density `ρ[i] · q⁺[i] · q⁻[i]` (probability per unit
volume). For the integrated reactive probability use
[`probability_reactive`](@ref).
"""
reactive_density(r::ReactiveTransition) = r.ρ .* r.qplus .* r.qminus

"""
$(TYPEDSIGNATURES)

Total reactive probability `∫ ρ · q⁺ · q⁻ dV` integrated over the grid.
"""
function probability_reactive(r::ReactiveTransition)::Float64
    v = cell_volume(r.generator.grid)
    s = 0.0
    @inbounds for i in eachindex(r.ρ)
        s += r.ρ[i] * r.qplus[i] * r.qminus[i]
    end
    return v * s
end

"""
$(TYPEDSIGNATURES)

Total probability that the most recent metastable visit was set A:
`∫ ρ · q⁻ dV`.
"""
function probability_last_A(r::ReactiveTransition)::Float64
    return cell_volume(r.generator.grid) * dot(r.ρ, r.qminus)
end

"""
$(TYPEDSIGNATURES)

Reactive probability current as `(J_nodes, J_faces)`:

- `J_faces :: NTuple{D, Array{Float64, D}}` — net signed flux across each
  axis-aligned face, divided by the corresponding axis spacing. Element
  `k` has shape `nbox` along all axes except axis `k`, where it has
  `nbox[k] - 1` entries for **reflecting / absorbing** boundary conditions
  (one per interior face), or `nbox[k]` entries for **periodic** boundary
  conditions (the extra slot is the wraparound face from cell `nbox[k]`
  back to cell `1`).
- `J_nodes :: NTuple{D, Array{Float64, D}}` — face fluxes averaged onto
  cell centers along each axis; element `k` has shape `nbox`.

The face flux from `n` to `m` is `μ[n] · q⁻[n] · Q[n, m] · q⁺[m]`, where
`μ = ρ · v` and `Q` is the rate matrix.

For non-equilibrium systems the current splits into a reversible
(gradient-flow) part and an irreversible (cyclic) part — see
[`reactive_current_reversible`](@ref) and
[`reactive_current_irreversible`](@ref). The three currents satisfy
`J_full = J_rev + J_irr` componentwise.
"""
function reactive_current(r::ReactiveTransition)
    return _reactive_current_from_Q(
        r.generator.Q, r.ρ, r.qplus, r.qminus, r.generator.grid, r.generator.bc,
    )
end

"""
$(TYPEDSIGNATURES)

Reversible (gradient-flow) part of the reactive current. Built from the
symmetric part of the generator, `Q_sym = (Q + Q̃)/2`, where `Q̃` is the
discrete adjoint with respect to the invariant density. For a reversible
system this equals the full [`reactive_current`](@ref).

Returns `(J_nodes, J_faces)` with the same shapes as
[`reactive_current`](@ref).
"""
function reactive_current_reversible(r::ReactiveTransition)
    Qadj = _adjoint_generator(r.generator.Q, r.ρ)
    Qsym = (r.generator.Q + Qadj) ./ 2
    return _reactive_current_from_Q(
        Qsym, r.ρ, r.qplus, r.qminus, r.generator.grid, r.generator.bc,
    )
end

"""
$(TYPEDSIGNATURES)

Irreversible (cyclic / time-asymmetric) part of the reactive current.
Built from the antisymmetric part of the generator,
`Q_antisym = (Q − Q̃)/2`. Vanishes for systems satisfying detailed
balance; for non-equilibrium systems it exposes the circulating component
of the reactive flow.

Returns `(J_nodes, J_faces)` with the same shapes as
[`reactive_current`](@ref).
"""
function reactive_current_irreversible(r::ReactiveTransition)
    Qadj = _adjoint_generator(r.generator.Q, r.ρ)
    Qanti = (r.generator.Q - Qadj) ./ 2
    return _reactive_current_from_Q(
        Qanti, r.ρ, r.qplus, r.qminus, r.generator.grid, r.generator.bc,
    )
end

function _reactive_current_from_Q(
        Q::SparseMatrixCSC{Float64, Int},
        ρ::Vector{Float64},
        qplus::Vector{Float64},
        qminus::Vector{Float64},
        grid::CartesianGrid{D},
        bc::Tuple,
    ) where {D}
    nbox = grid.nbox
    LI = LinearIndices(nbox)
    v = cell_volume(grid)

    # Periodic axes get one extra face slot (the wraparound face).
    J_faces = ntuple(D) do k
        len_k = bc[k] isa Periodic ? nbox[k] : nbox[k] - 1
        zeros(Float64, ntuple(d -> d == k ? len_k : nbox[d], D)...)
    end

    @inbounds for k in 1:D
        hk = grid.h[k]
        Jk = J_faces[k]
        periodic_k = bc[k] isa Periodic
        for I in CartesianIndices(nbox)
            if I[k] == nbox[k] && !periodic_k
                continue
            end
            j_k = I[k] == nbox[k] ? 1 : I[k] + 1            # wraparound index
            J = CartesianIndex(ntuple(d -> d == k ? j_k : I[d], D))
            n = LI[I]
            m = LI[J]
            μn = ρ[n] * v
            μm = ρ[m] * v
            flux_nm = μn * qminus[n] * Q[n, m] * qplus[m]
            flux_mn = μm * qminus[m] * Q[m, n] * qplus[n]
            Jk[I] = (flux_nm - flux_mn) / hk
        end
    end

    J_nodes = ntuple(_ -> zeros(Float64, nbox...), D)
    @inbounds for k in 1:D
        Jface = J_faces[k]
        Jnode = J_nodes[k]
        nk = nbox[k]
        periodic_k = bc[k] isa Periodic
        for I in CartesianIndices(nbox)
            s = 0.0
            c = 0
            # Previous face: between cell I-e_k and cell I.
            if I[k] > 1
                Iprev = CartesianIndex(ntuple(d -> d == k ? I[d] - 1 : I[d], D))
                s += Jface[Iprev]
                c += 1
            elseif periodic_k
                # Wraparound: prev face stored at index nbox[k] of axis k.
                Iprev = CartesianIndex(ntuple(d -> d == k ? nk : I[d], D))
                s += Jface[Iprev]
                c += 1
            end
            # Next face: between cell I and cell I+e_k.
            if I[k] < nk
                s += Jface[I]
                c += 1
            elseif periodic_k
                # Wraparound: next face stored at index nbox[k] (== I[k]).
                s += Jface[I]
                c += 1
            end
            c > 0 && (Jnode[I] = s / c)
        end
    end

    return J_nodes, J_faces
end
