"""
$(TYPEDEF)

A `D`-dimensional uniform axis-aligned Cartesian grid with coordinates
stored in the floating-point type `T` (default `Float64`).

# Fields
$(TYPEDFIELDS)

# Construction
```julia
CartesianGrid((u_lo, u_hi, N_u), (v_lo, v_hi, N_v), ...)            # Float64
CartesianGrid{Double64}((u_lo, u_hi, N_u), (v_lo, v_hi, N_v), ...)  # extended
```
Each axis is `(lo, hi, N)` and creates `N` equal-width cells. The
storage type `T` propagates into [`DiffusionGenerator`](@ref) and the
resulting rate matrix.
"""
struct CartesianGrid{D, T <: AbstractFloat}
    "Number of cells along each axis."
    nbox::NTuple{D, Int}
    "Cell width along each axis."
    h::NTuple{D, T}
    "Cell-center coordinates along each axis (length `N_k`)."
    centers::NTuple{D, LinRange{T, Int}}
    "Cell-edge coordinates along each axis (length `N_k + 1`)."
    edges::NTuple{D, LinRange{T, Int}}
end

function CartesianGrid{T}(
        axes::Vararg{Tuple{Real, Real, Integer}, D}
    ) where {D, T <: AbstractFloat}
    D >= 1 || throw(ArgumentError("CartesianGrid needs at least one axis"))
    nbox = ntuple(k -> Int(axes[k][3]), D)
    @inbounds for k in 1:D
        nbox[k] >= 2 || throw(ArgumentError("axis $k needs ≥ 2 cells, got $(nbox[k])"))
        axes[k][2] > axes[k][1] ||
            throw(ArgumentError("axis $k: hi=$(axes[k][2]) must exceed lo=$(axes[k][1])"))
    end
    h = ntuple(k -> (T(axes[k][2]) - T(axes[k][1])) / nbox[k], D)
    edges = ntuple(
        k -> LinRange(T(axes[k][1]), T(axes[k][2]), nbox[k] + 1), D
    )
    centers = ntuple(
        k -> LinRange(
            T(axes[k][1]) + h[k] / 2,
            T(axes[k][2]) - h[k] / 2,
            nbox[k],
        ),
        D,
    )
    return CartesianGrid{D, T}(nbox, h, centers, edges)
end

CartesianGrid(axes::Vararg{Tuple{Real, Real, Integer}, D}) where {D} =
    CartesianGrid{Float64}(axes...)

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

Center coordinates of cell `I` of `grid`.
"""
@inline function cell_center(
        grid::CartesianGrid{D, T}, I::CartesianIndex{D}
    ) where {D, T}
    return SVector{D, T}(ntuple(d -> grid.centers[d][I[d]], D))
end
