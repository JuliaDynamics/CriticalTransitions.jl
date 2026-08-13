"""
$(TYPEDEF)

Sub-cell back-reference used to reconstruct instantons from a
[`QuasiPotential`](@ref) grid. For one-point (vertex) updates `v1 == v0` and
`s == NaN32`. For two-point (edge) updates the predecessor lies at
`(1 - s) * v0 + s * v1` with `s ∈ [0, 1]`.

# Fields
$(TYPEDFIELDS)
"""
struct BackRef{D}
    v0::CartesianIndex{D}
    v1::CartesianIndex{D}
    s::Float32
end

BackRef{D}() where {D} =
    BackRef{D}(zero(CartesianIndex{D}), zero(CartesianIndex{D}), NaN32)

"""
$(TYPEDEF)

Discrete quasipotential field `U_A(x)` computed by [`quasipotential`](@ref).
`U` is `+Inf` on cells the sweep did not reach. `back_pointer` encodes the
winning sub-cell predecessor per cell (zeroed for the source and unreached
cells).

# Fields
$(TYPEDFIELDS)
"""
struct QuasiPotential{D, T <: AbstractFloat}
    U::Array{T, D}
    back_pointer::Array{BackRef{D}, D}
    source::CartesianIndex{D}
    grid::CartesianGrid{D, T}
end

@enum Status::UInt8 _UNKNOWN _CONSIDERED _FRONT _ACCEPTED

struct _OLIMState{D, T, H}
    U::Array{T, D}
    back_pointer::Array{BackRef{D}, D}
    status::Array{Status, D}
    front::Array{Bool, D}
    heap::H
    handles::Vector{Int}
    nbox::NTuple{D, Int}
    # Per-cell drift cache for the additive-noise line integral: `b_c[I] = b(c_I)`,
    # `q_c[I] = |b(c_I)|²_Q`, `sq_c[I] = |b(c_I)|_Q` at each cell center `c_I`.
    # Populated once by `_fill_node_cache!` (additive noise only). For multiplicative
    # noise the metric is not constant, so these are allocated empty (zero-size,
    # never indexed) to avoid a full-grid allocation that would never be read.
    b_c::Array{SVector{D, T}, D}
    q_c::Array{T, D}
    sq_c::Array{T, D}
end

# `cache_drift` allocates the full-grid additive-noise drift cache; pass `false`
# for multiplicative noise to allocate empty placeholders (see `_fill_node_cache!`).
function _OLIMState(
        grid::CartesianGrid{D, T}, ::Type{T}, cache_drift::Bool = true,
    ) where {D, T}
    nbox = grid.nbox
    U = fill(T(Inf), nbox)
    bp = fill(BackRef{D}(), nbox)
    status = fill(_UNKNOWN, nbox)
    front = fill(false, nbox)
    heap = MutableBinaryHeap{Tuple{T, Int}, FasterForward}()
    # Pre-size the heap's internal storage to the maximum possible (one slot
    # per cell). The heap will never need to grow during the sweep.
    sizehint!(heap.nodes, prod(nbox))
    sizehint!(heap.node_map, prod(nbox))
    handles = zeros(Int, prod(nbox))
    cbox = cache_drift ? nbox : ntuple(_ -> 0, Val(D))
    b_c = fill(zero(SVector{D, T}), cbox)
    q_c = zeros(T, cbox)
    sq_c = zeros(T, cbox)
    return _OLIMState{D, T, typeof(heap)}(
        U, bp, status, front, heap, handles, nbox, b_c, q_c, sq_c,
    )
end

@inline _linear(x::CartesianIndex{D}, nbox::NTuple{D, Int}) where {D} =
    LinearIndices(nbox)[x]
