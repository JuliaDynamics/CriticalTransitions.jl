"""
$(TYPEDEF)

A `D`-dimensional uniform axis-aligned Cartesian grid. Fully concrete and
parameterised over the dimension `D`, so all subsequent operations (chain
assembly, Cartesian indexing, accessors) specialise at compile time.

# Fields
$(TYPEDFIELDS)

# Construction
```julia
CartesianGrid((u_lo, u_hi, N_u), (v_lo, v_hi, N_v), ...)
```
Each axis is given as `(lo, hi, N)` and produces `N` cells of equal width
covering `[lo, hi]`.
"""
struct CartesianGrid{D}
    "Number of cells along each axis."
    nbox::NTuple{D, Int}
    "Cell width along each axis."
    h::NTuple{D, Float64}
    "Cell-center coordinates along each axis (length `N_k`)."
    centers::NTuple{D, LinRange{Float64, Int}}
    "Cell-edge coordinates along each axis (length `N_k + 1`)."
    edges::NTuple{D, LinRange{Float64, Int}}
end

function CartesianGrid(axes::Vararg{Tuple{Real, Real, Integer}, D}) where {D}
    D >= 1 || throw(ArgumentError("CartesianGrid needs at least one axis"))
    nbox = ntuple(k -> Int(axes[k][3]), D)
    @inbounds for k in 1:D
        nbox[k] >= 2 || throw(ArgumentError("axis $k needs ≥ 2 cells, got $(nbox[k])"))
        axes[k][2] > axes[k][1] ||
            throw(ArgumentError("axis $k: hi=$(axes[k][2]) must exceed lo=$(axes[k][1])"))
    end
    h = ntuple(k -> (Float64(axes[k][2]) - Float64(axes[k][1])) / nbox[k], D)
    edges = ntuple(
        k -> LinRange(Float64(axes[k][1]), Float64(axes[k][2]), nbox[k] + 1), D
    )
    centers = ntuple(
        k -> LinRange(
            Float64(axes[k][1]) + h[k] / 2,
            Float64(axes[k][2]) - h[k] / 2,
            nbox[k],
        ),
        D,
    )
    return CartesianGrid{D}(nbox, h, centers, edges)
end

"""
$(TYPEDSIGNATURES)

Total number of cells in `grid`.
"""
@inline ncells(grid::CartesianGrid) = prod(grid.nbox)

"""
$(TYPEDSIGNATURES)

Volume of one cell (product of axis spacings).
"""
@inline cell_volume(grid::CartesianGrid) = prod(grid.h)

"""
$(TYPEDSIGNATURES)

`SVector{D, Float64}` of the center coordinates of cell `I`.
"""
@inline function cell_center(grid::CartesianGrid{D}, I::CartesianIndex{D}) where {D}
    return SVector{D, Float64}(ntuple(d -> grid.centers[d][I[d]], D))
end
