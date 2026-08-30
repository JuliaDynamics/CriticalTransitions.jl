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
    return _shift_in_bounds(x, δ, nbox, ntuple(_ -> false, Val(D)))
end

@inline function _shift_in_bounds(
        x::CartesianIndex{D}, δ::CartesianIndex{D}, nbox::NTuple{D, Int},
        periodic::NTuple{D, Bool},
    ) where {D}
    idx = ntuple(Val(D)) do d
        i = x[d] + δ[d]
        if periodic[d]
            mod1(i, nbox[d])
        else
            (i < 1 || i > nbox[d]) && return 0
            i
        end
    end
    @inbounds for d in 1:D
        idx[d] == 0 && return nothing
    end
    return CartesianIndex(idx)
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
    return _is_chebyshev_adjacent(
        p, q, ntuple(_ -> typemax(Int), Val(D)), ntuple(_ -> false, Val(D)),
    )
end

@inline function _wrapped_index_distance(a::Int, b::Int, n::Int, periodic::Bool)
    δ = abs(a - b)
    return periodic ? min(δ, n - δ) : δ
end

@inline function _is_chebyshev_adjacent(
        p::CartesianIndex{D}, q::CartesianIndex{D},
        nbox::NTuple{D, Int}, periodic::NTuple{D, Bool},
    ) where {D}
    p == q && return false
    @inbounds for d in 1:D
        _wrapped_index_distance(p[d], q[d], nbox[d], periodic[d]) > 1 &&
            return false
    end
    return true
end

@inline function _in_circular_stencil(
        x::CartesianIndex{D}, y::CartesianIndex{D}, ::Val{K}
    ) where {D, K}
    return _in_circular_stencil(
        x, y, Val(K), ntuple(_ -> typemax(Int), Val(D)),
        ntuple(_ -> false, Val(D)),
    )
end

@inline function _in_circular_stencil(
        x::CartesianIndex{D}, y::CartesianIndex{D}, ::Val{K},
        nbox::NTuple{D, Int}, periodic::NTuple{D, Bool},
    ) where {D, K}
    s = 0
    @inbounds for d in 1:D
        δ = _wrapped_index_distance(y[d], x[d], nbox[d], periodic[d])
        s += δ * δ
    end
    return s != 0 && s <= K * K
end

"""
Return the center of `I`, lifted by whole periods so it is the nearest periodic
image to the center of `anchor`. Nonperiodic coordinates are unchanged.
"""
@inline function _cell_center_near(
        grid::CartesianGrid{D, T}, I::CartesianIndex{D},
        anchor::CartesianIndex{D}, periodic::NTuple{D, Bool},
    ) where {D, T}
    c = cell_center(grid, I)
    a = cell_center(grid, anchor)
    return SVector{D, T}(
        ntuple(Val(D)) do d
            periodic[d] || return c[d]
            period = grid.h[d] * grid.nbox[d]
            c[d] + round((a[d] - c[d]) / period) * period
        end
    )
end
