"""
$(TYPEDEF)

Transition-path-theory result for transitions from set `A` to set `B`.

`ReactiveTransition` stores the generator, set masks, invariant density,
forward and backward committors, discrete adjoint generator, and reactive
rate. It represents the stationary ensemble of trajectories that most
recently left `A` and will hit `B` before returning to `A`.

# Fields
$(TYPEDFIELDS)

# Construction
```julia
ReactiveTransition(gen::DiffusionGenerator, A, B; alg=UMFPACKFactorization(), reverse=nothing)
ReactiveTransition(sys::CoupledSDEs, grid, A, B; bc=Reflecting(), reverse=nothing)
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


```julia
res = ReactiveTransition(gen, x -> x[1] < -1, x -> x[1] > 1)
k = reactive_rate(res)
J_nodes, J_faces = reactive_current(res)
```

Default `alg` is `UMFPACKFactorization()` (sparse direct LU) for all
internal solves. Pass any `SciMLLinearSolveAlgorithm` (e.g.
`KrylovJL_GMRES()`) for an iterative solver — useful for large grids
where direct LU runs out of memory. Further kwargs flow to
`LinearSolve.solve`.

See also [`forward_committor`](@ref), [`backward_committor`](@ref),
[`reactive_density`](@ref), and [`reactive_current`](@ref).
"""
struct ReactiveTransition{D, BC, T <: AbstractFloat}
    "The diffusion generator the analysis was built on."
    generator::DiffusionGenerator{D, BC, T}
    "Boolean mask selecting cells in set A (length `ncells(grid)`)."
    A_mask::BitVector
    "Boolean mask selecting cells in set B (length `ncells(grid)`)."
    B_mask::BitVector
    "Invariant probability density (`dot(ρ, weights) = 1`)."
    ρ::Vector{T}
    "Forward committor `q⁺` (probability of reaching B before A)."
    qplus::Vector{T}
    "Backward committor `q⁻`."
    qminus::Vector{T}
    "Discrete adjoint generator w.r.t. `ρ`; used by `reactive_current_{reversible,irreversible}`."
    Qadj::SparseMatrixCSC{T, Int}
    "A → B reactive transition rate (cached)."
    rate::T
    "True if `qminus` came from a user-supplied physical reverse generator; false if from the discrete adjoint."
    physical_reverse::Bool
end

function ReactiveTransition(
        gen::DiffusionGenerator{D, BC, T},
        A,
        B;
        alg = UMFPACKFactorization(),
        reverse::Union{Nothing, DiffusionGenerator{D}} = nothing,
        kwargs...,
    ) where {D, BC, T}
    any(b -> b isa Absorbing, gen.bc) && throw(
        ArgumentError(
            "ReactiveTransition is undefined for absorbing boundary " *
                "conditions: the chain has mass loss and admits no nontrivial " *
                "invariant density. Use Reflecting() or Periodic() instead.",
        ),
    )
    grid = gen.grid
    A_mask, B_mask = _committor_masks(grid, A, B)

    Q = gen.Q
    weights = fill(cell_volume(grid), ncells(gen))

    ρ = stationary_distribution(gen, alg; kwargs...)
    qplus = _forward_committor(Q, A_mask, B_mask, alg; kwargs...)
    Qadj = _adjoint_generator(Q, ρ)

    if reverse === nothing
        qminus = _backward_committor_explicit(Qadj, A_mask, B_mask, alg; kwargs...)
        physical_reverse = false
    else
        _check_reverse_generator(gen, reverse)
        qminus = _backward_committor_explicit(reverse.Q, A_mask, B_mask, alg; kwargs...)
        physical_reverse = true
    end

    rate = _reactive_rate(Q, ρ, qplus, qminus, weights, A_mask)
    return ReactiveTransition{D, typeof(gen.bc), T}(
        gen, A_mask, B_mask, ρ, qplus, qminus, Qadj, rate, physical_reverse,
    )
end

function ReactiveTransition(
        sys::CoupledSDEs,
        grid::CartesianGrid{D, T},
        A,
        B;
        alg = UMFPACKFactorization(),
        reverse::Union{Nothing, CoupledSDEs} = nothing,
        bc = Reflecting(),
        kwargs...,
    ) where {D, T}
    gen = DiffusionGenerator(sys, grid; bc = bc)
    gen_rev =
        reverse === nothing ? nothing : DiffusionGenerator(reverse, grid; bc = bc)
    return ReactiveTransition(gen, A, B; alg, reverse = gen_rev, kwargs...)
end

@inline forward_committor(r::ReactiveTransition) = r.qplus
@inline backward_committor(r::ReactiveTransition) = r.qminus
@inline stationary_distribution(r::ReactiveTransition) = r.ρ

"""
$(TYPEDSIGNATURES)

A → B **reactive transition rate** ``k_{AB}``: events per unit time, the
rate at which probability leaves `A` along reactive trajectories.

In CTMC form,

    k_{AB} = ∑_{i ∈ A, j ∉ A} ρ[i] · v · q⁻[i] · Q[i, j] · q⁺[j] ,

with `v = cell_volume(grid)`.
"""
@inline reactive_rate(r::ReactiveTransition) = r.rate

"""
$(TYPEDSIGNATURES)

Per-cell reactive density `ρ[i] * q⁺[i] * q⁻[i]`. This is a density per unit volume.
For the integrated reactive probability use
[`probability_reactive`](@ref).
"""
reactive_density(r::ReactiveTransition) = r.ρ .* r.qplus .* r.qminus

"""
$(TYPEDSIGNATURES)

Total reactive probability integrated over the grid.
This is the stationary probability that the process is currently on a
reactive segment from `A` to `B`.

See also [`reactive_density`](@ref).
"""
function probability_reactive(r::ReactiveTransition{D, BC, T})::T where {D, BC, T}
    v = cell_volume(r.generator.grid)
    s = zero(T)
    @inbounds for i in eachindex(r.ρ)
        s += r.ρ[i] * r.qplus[i] * r.qminus[i]
    end
    return v * s
end

"""
$(TYPEDSIGNATURES)

Total probability that the most recent metastable visit was set `A`.

The complementary probability is the chance that the most recent visit was
`B`.
"""
function probability_last_A(r::ReactiveTransition{D, BC, T})::T where {D, BC, T}
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
[`reactive_current_irreversible`](@ref) to split gradient-flow and cyclic
components.
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
"""
function reactive_current_reversible(r::ReactiveTransition)
    return _reactive_current_from_Qs(
        (r.generator.Q, r.Qadj), (0.5, 0.5),
        r.ρ, r.qplus, r.qminus, r.generator.grid, r.generator.bc,
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
    return _reactive_current_from_Qs(
        (r.generator.Q, r.Qadj), (0.5, -0.5),
        r.ρ, r.qplus, r.qminus, r.generator.grid, r.generator.bc,
    )
end

function _reactive_current_from_Q(
        Q::SparseMatrixCSC{T, Int},
        ρ::Vector{T},
        qplus::Vector{T},
        qminus::Vector{T},
        grid::CartesianGrid{D, T},
        bc::Tuple,
    ) where {D, T}
    return _reactive_current_from_Qs((Q,), (one(T),), ρ, qplus, qminus, grid, bc)
end

@inline function _edge_rate(Qs::Tuple, coeffs::Tuple, i::Int, j::Int)
    T = eltype(Qs[1])
    rate = zero(T)
    @inbounds for k in eachindex(Qs)
        rate += coeffs[k] * Qs[k][i, j]
    end
    return rate
end

function _reactive_current_from_Qs(
        Qs::Tuple,
        coeffs::Tuple,
        ρ::Vector{T},
        qplus::Vector{T},
        qminus::Vector{T},
        grid::CartesianGrid{D, T},
        bc::Tuple,
    ) where {D, T}
    nbox = grid.nbox
    LI = LinearIndices(nbox)
    v = cell_volume(grid)

    # Periodic axes get one extra face slot (the wraparound face).
    J_faces = ntuple(D) do k
        len_k = bc[k] isa Periodic ? nbox[k] : nbox[k] - 1
        zeros(T, ntuple(d -> d == k ? len_k : nbox[d], D)...)
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
            flux_nm = μn * qminus[n] * _edge_rate(Qs, coeffs, n, m) * qplus[m]
            flux_mn = μm * qminus[m] * _edge_rate(Qs, coeffs, m, n) * qplus[n]
            Jk[I] = (flux_nm - flux_mn) / hk
        end
    end

    J_nodes = ntuple(_ -> zeros(T, nbox...), D)
    @inbounds for k in 1:D
        Jface = J_faces[k]
        Jnode = J_nodes[k]
        nk = nbox[k]
        periodic_k = bc[k] isa Periodic
        for I in CartesianIndices(nbox)
            s = zero(T)
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

"""
$(TYPEDSIGNATURES)

A → B reactive transition rate from a CTMC generator `G` (positive
off-diagonals = transition rates):

    k_{AB} = ∑_{i ∈ A, j ∉ A} ρ[i] · v[i] · q⁻[i] · G[i, j] · q⁺[j] .

Computed in `O(nnz(G))` by walking sparse columns once.
"""
function _reactive_rate(
        G::SparseMatrixCSC{T, Int},
        ρ::Vector{T},
        qplus::Vector{T},
        qminus::Vector{T},
        weights::Vector{T},
        A_mask::BitVector,
    )::T where {T <: AbstractFloat}
    rv = rowvals(G)
    nz = nonzeros(G)
    rate = zero(T)
    @inbounds for col in 1:size(G, 2)
        A_mask[col] && continue
        qpc = qplus[col]
        for p in nzrange(G, col)
            row = rv[p]
            row != col || continue
            A_mask[row] || continue
            rate += ρ[row] * weights[row] * qminus[row] * nz[p] * qpc
        end
    end
    return rate
end
