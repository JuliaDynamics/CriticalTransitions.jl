# Circular-stencil offsets, sorted by squared radius. Returned as a `const` `Vector`
# (built once per `(K, D)`): iterating the cached vector is linear, whereas the old
# tuple form hit Julia's large-tuple iteration cliff for the wide stencils that
# anisotropic problems require (e.g. K=12 -> a ~450-element tuple). The same vector
# instance is shared across all calls with the same `(K, D)`; callers iterate it
# read-only and must never mutate it.
@generated function _stencil_offsets(::Val{K}, ::Val{D}) where {K, D}
    offsets = CartesianIndex{D}[]
    box = ntuple(_ -> -K:K, D)
    for I in Iterators.product(box...)
        all(==(0), I) && continue
        sum(abs2, I) <= K * K || continue
        push!(offsets, CartesianIndex(I))
    end
    sort!(offsets; by = t -> sum(abs2, Tuple(t)))
    return :($offsets)
end

@inline function _shift_in_bounds(
        x::CartesianIndex{D}, δ::CartesianIndex{D}, nbox::NTuple{D, Int}
    ) where {D}
    @inbounds for d in 1:D
        i = x[d] + δ[d]
        (i < 1 || i > nbox[d]) && return nothing
    end
    return x + δ
end

@generated function _chebyshev_neighbors(::Val{D}) where {D}
    offsets = NTuple{D, Int}[]
    box = ntuple(_ -> -1:1, D)
    for I in Iterators.product(box...)
        all(==(0), I) && continue
        push!(offsets, I)
    end
    cis = [:(CartesianIndex($(t...))) for t in offsets]
    return Expr(:tuple, cis...)
end

@inline function _is_chebyshev_adjacent(
        p::CartesianIndex{D}, q::CartesianIndex{D}
    ) where {D}
    p == q && return false
    @inbounds for d in 1:D
        abs(p[d] - q[d]) > 1 && return false
    end
    return true
end

@inline function _in_circular_stencil(
        x::CartesianIndex{D}, y::CartesianIndex{D}, ::Val{K}
    ) where {D, K}
    s = 0
    @inbounds for d in 1:D
        δ = y[d] - x[d]
        s += δ * δ
    end
    return s != 0 && s <= K * K
end
