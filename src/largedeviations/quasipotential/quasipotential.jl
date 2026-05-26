"""
    default_K(grid::CartesianGrid) -> Int

Heuristic accepted-band radius (in grid cells):
`max(5, round(Int, sqrt(minimum(grid.nbox))))`, capped at `32`.
"""
@inline default_K(grid::CartesianGrid) =
    clamp(round(Int, sqrt(minimum(grid.nbox))), 5, 32)

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
* `verbose::Bool      = false`
* `show_progress::Bool = true`

See also: [`QuasiPotential`](@ref), [`BackRef`](@ref).
"""
function quasipotential(
        sys::CoupledSDEs{IIP, D, I, P},
        grid::CartesianGrid{D, T},
        attractor::AbstractVector{<:Real};
        band_radius::Int    = default_K(grid),
        near_source_layers::Int = 3,
        verbose::Bool       = false,
        show_progress::Bool = true,
    ) where {IIP, D, I, P, T}
    length(attractor) == D || throw(
        DimensionMismatch(
            "attractor has length $(length(attractor)) but sys has D=$D",
        ),
    )
    proper_FW_system(sys)
    D > 4 && @warn "quasipotential in D=$D: per-axis grid resolution will be coarse" maxlog=1
    x_A = SVector{D, T}(attractor)
    return _quasipotential_impl(
        sys, grid, x_A,
        Val(band_radius), Val(near_source_layers),
        verbose, show_progress,
    )
end

function _quasipotential_impl(
        sys::CoupledSDEs{IIP, D, I, P},
        grid::CartesianGrid{D, T},
        x_A::SVector{D, T},
        ::Val{K}, ::Val{K_seed},
        verbose::Bool, show_progress::Bool,
    ) where {IIP, D, I, P, T, K, K_seed}
    src = _source_cell(grid, x_A)
    state = _OLIMState(grid, T)
    L = _geometric_lagrangian(sys, T)
    _sweep!(
        state, grid, src, sys, L, Val(K), Val(K_seed);
        verbose = verbose, show_progress = show_progress,
    )
    return QuasiPotential{D, T}(state.U, state.back_pointer, src, grid)
end
