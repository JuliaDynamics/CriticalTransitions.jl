"""
    default_K(grid::CartesianGrid{D}) -> Int

Heuristic accepted-band radius in grid cells, dimension-aware. Per-cell sweep
work scales as `(2K+1)^D`, so the floor and cap depend on `D`:

| `D`   | floor | cap |
|-------|-------|-----|
| 1, 2  |   5   |  32 |
| 3     |   2   |   8 |
| ≥ 4   |   2   |   4 |

Within these bounds the radius is `round(Int, sqrt(minimum(grid.nbox)))`.
"""
@inline function default_K(grid::CartesianGrid{D}) where {D}
    floor_K, cap_K = D <= 2 ? (5, 32) : D == 3 ? (2, 8) : (2, 4)
    return clamp(round(Int, sqrt(minimum(grid.nbox))), floor_K, cap_K)
end

"""
    default_regularization(grid::CartesianGrid) -> T

Default regularization for rank-1 (single noiseless coordinate) diffusion: the amount
added to the noiseless diagonal of the trace-normalized diffusion so the metric is
invertible. Proportional to the noiseless cell fraction (`1 / minimum(grid.nbox)`) so it
vanishes under refinement, with the constant calibrated to about `0.04` at `80` cells
across. Ignored for full-rank systems; override with the `regularization` keyword of
[`quasipotential`](@ref).
"""
@inline default_regularization(grid::CartesianGrid{D, T}) where {D, T} =
    T(3.2) / minimum(grid.nbox)

function _source_cell(
        grid::CartesianGrid{D, T}, attractor::SVector{D, T}
    ) where {D, T}
    idx = ntuple(D) do d
        c = grid.centers[d]
        lo = c[1] - grid.h[d] / 2
        hi = c[end] + grid.h[d] / 2
        (lo <= attractor[d] <= hi) ||
            throw(
            ArgumentError(
                "attractor[$d] = $(attractor[d]) outside grid axis $d ($lo, $hi)",
            ),
        )
        i = round(Int, (attractor[d] - c[1]) / grid.h[d]) + 1
        clamp(i, 1, length(c))
    end
    return CartesianIndex(idx)
end

"""
    quasipotential(sys, grid, attractor; band_radius, near_source_layers,
                   verbose, show_progress) -> QuasiPotential{D}

Compute the Freidlin-Wentzell quasipotential field `U_A(x)` from `attractor`
using the Ordered Line Integral Method (Dahiya and Cameron 2018). The state
dimension `D` is taken from `sys::CoupledSDEs{IIP, D, I, P}` and must match
`grid::CartesianGrid{D}`. A warning is emitted for `D > 4`.

# Keyword arguments
* `band_radius::Int  = default_K(grid)`: accepted-band radius in grid cells.
* `near_source_layers::Int = 3`: size of the analytic CARE seed box;
  `0` disables analytic seeding.
* `regularization::Real = default_regularization(grid)`: for systems with one
  noiseless coordinate (rank-1 diffusion), the amount added to that coordinate's
  diagonal of the trace-normalized diffusion to make the metric invertible. Smaller
  values are more accurate but stiffer; reduce it as the grid is refined. Ignored for
  full-rank systems.
* `verbose::Bool      = false`
* `show_progress::Bool = true`

See also: [`QuasiPotential`](@ref), [`BackRef`](@ref).
"""
function quasipotential(
        sys::CoupledSDEs{IIP, D, I, P},
        grid::CartesianGrid{D, T},
        attractor::AbstractVector{<:Real};
        band_radius::Int = default_K(grid),
        near_source_layers::Int = 3,
        regularization::Real = default_regularization(grid),
        verbose::Bool = false,
        show_progress::Bool = true,
    ) where {IIP, D, I, P, T}
    length(attractor) == D || throw(
        DimensionMismatch(
            "attractor has length $(length(attractor)) but sys has D=$D",
        ),
    )
    proper_FW_system(sys)
    D > 4 && @warn "quasipotential in D=$D: per-axis grid resolution will be coarse" maxlog = 1
    x_A = SVector{D, T}(attractor)
    return _quasipotential_impl(
        sys, grid, x_A,
        Val(band_radius), Val(near_source_layers),
        T(regularization), verbose, show_progress,
    )
end

function _quasipotential_impl(
        sys::CoupledSDEs{IIP, D, I, P},
        grid::CartesianGrid{D, T},
        x_A::SVector{D, T},
        ::Val{K}, ::Val{K_seed},
        regularization::T, verbose::Bool, show_progress::Bool,
    ) where {IIP, D, I, P, T, K, K_seed}
    src = _source_cell(grid, x_A)
    state = _OLIMState(grid, T)
    L = _geometric_lagrangian(sys, T; regularization = regularization)
    _sweep!(
        state, grid, src, sys, L, Val(K), Val(K_seed);
        verbose = verbose, show_progress = show_progress,
    )
    return QuasiPotential{D, T}(state.U, state.back_pointer, src, grid)
end
