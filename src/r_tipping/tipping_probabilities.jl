"""
    tipping_probabilities(
        BoA_before::ArrayBasinsOfAttraction,
        BoA_after::ArrayBasinsOfAttraction,
    ) → P
    tipping_probabilities(basins_before, basins_after) → P

Return the basin-overlap tipping probabilities before and after a parameter change or
time-dependent forcing, following [Kaszas2019](@cite).

The basin arrays must have identical sizes and contain integer attractor labels. Let
``\mathcal{B}_i(p)`` denote the basin of attraction of attractor ``A_i`` at parameters
``p``. For a change ``p_- \to p_+``, the matrix entry is

```math
P(A_i \to A_j \mid p_- \to p_+) =
\frac{|\mathcal{B}_j(p_+) \cap \mathcal{B}_i(p_-)|}{|\mathcal{B}_i(p_-)|}.
```

Rows are ordered by the sorted unique labels in `basins_before`, and columns by the sorted
unique labels in `basins_after`. If the label `-1` occurs (trajectories that diverge), it is
placed last. For the standard consecutive labels `1:n`, `P[i, j]` is therefore the tipping
probability from attractor `i` to attractor `j`.

`ArrayBasinsOfAttraction` inputs are supported for convenience; their `.basins` arrays are
used directly and must describe the same grid points in the same order.
"""
function tipping_probabilities(
    basins_before::AbstractArray{<:Integer}, basins_after::AbstractArray{<:Integer}
)
    size(basins_before) == size(basins_after) ||
        throw(DimensionMismatch("basin arrays must have identical sizes"))

    before_ids = _tipping_probability_ids(basins_before)
    after_ids = _tipping_probability_ids(basins_after)
    before_index = Dict(id => i for (i, id) in pairs(before_ids))
    after_index = Dict(id => i for (i, id) in pairs(after_ids))

    counts = zeros(Int, length(before_ids), length(after_ids))
    totals = zeros(Int, length(before_ids))
    for (before, after) in zip(basins_before, basins_after)
        i = before_index[before]
        j = after_index[after]
        counts[i, j] += 1
        totals[i] += 1
    end

    P = zeros(Float64, size(counts))
    for i in axes(counts, 1), j in axes(counts, 2)
        P[i, j] = counts[i, j] / totals[i]
    end
    return P
end

function tipping_probabilities(
    BoA_before::ArrayBasinsOfAttraction, BoA_after::ArrayBasinsOfAttraction
)
    return tipping_probabilities(BoA_before.basins, BoA_after.basins)
end

function _tipping_probability_ids(basins)
    ids = sort!(collect(unique(basins)))
    divergent = findfirst(isequal(-1), ids)
    isnothing(divergent) || push!(ids, popat!(ids, divergent))
    return ids
end
