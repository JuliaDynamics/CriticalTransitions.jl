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
    reverse::Union{Nothing,DiffusionGenerator{D}}=nothing,
) where {D}
    grid = gen.grid
    A_mask = _to_mask(A, grid)
    B_mask = _to_mask(B, grid)
    any(A_mask) || throw(ArgumentError("set A is empty"))
    any(B_mask) || throw(ArgumentError("set B is empty"))
    any(A_mask .& B_mask) && throw(ArgumentError("sets A and B overlap"))

    Q = gen.Q
    weights = fill(cell_volume(grid), ncells(gen))

    ρ = _invariant_density(Q, weights)
    qplus = _forward_committor(Q, A_mask, B_mask)

    if reverse === nothing
        qminus = _backward_committor_adjoint(Q, ρ, A_mask, B_mask)
        physical_reverse = false
    else
        qminus = _backward_committor_explicit(reverse.Q, A_mask, B_mask)
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
    reverse::Union{Nothing,CoupledSDEs}=nothing,
) where {D}
    gen = DiffusionGenerator(sys, grid)
    gen_rev = reverse === nothing ? nothing : DiffusionGenerator(reverse, grid)
    return ReactiveTransition(gen, A, B; reverse=gen_rev)
end

# =====================================================================
# Cached accessors (constant-time field reads)
# =====================================================================

@inline forward_committor(r::ReactiveTransition) = r.qplus
@inline backward_committor(r::ReactiveTransition) = r.qminus
@inline stationary_distribution(r::ReactiveTransition) = r.ρ
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
  axis-aligned face, divided by the corresponding axis spacing. Element `k`
  has shape `nbox` with `nbox[k] - 1` along axis `k` (one entry per face).
- `J_nodes :: NTuple{D, Array{Float64, D}}` — face fluxes averaged onto
  cell centers along each axis; element `k` has shape `nbox`.

The face flux from `n` to `m` is `μ[n] · q⁻[n] · Q[n, m] · q⁺[m]`, where
`μ = ρ · v` and `Q` is the rate matrix.
"""
function reactive_current(r::ReactiveTransition{D}) where {D}
    grid = r.generator.grid
    nbox = grid.nbox
    LI = LinearIndices(nbox)
    v = cell_volume(grid)
    Q = r.generator.Q
    qplus = r.qplus
    qminus = r.qminus
    ρ = r.ρ

    J_faces = ntuple(
        k -> zeros(Float64, ntuple(d -> d == k ? nbox[d] - 1 : nbox[d], D)...), D
    )

    @inbounds for k in 1:D
        hk = grid.h[k]
        Jk = J_faces[k]
        for I in CartesianIndices(nbox)
            I[k] == nbox[k] && continue
            J = CartesianIndex(ntuple(d -> d == k ? I[d] + 1 : I[d], D))
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
        for I in CartesianIndices(nbox)
            s = 0.0
            c = 0
            if I[k] > 1
                Iprev = CartesianIndex(ntuple(d -> d == k ? I[d] - 1 : I[d], D))
                s += Jface[Iprev]
                c += 1
            end
            if I[k] < nk
                s += Jface[I]
                c += 1
            end
            c > 0 && (Jnode[I] = s / c)
        end
    end

    return J_nodes, J_faces
end
