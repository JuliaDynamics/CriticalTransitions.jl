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

function _source_cells(
        grid::CartesianGrid{D, T},
        attractor::AbstractVector{<:AbstractVector},
    ) where {D, T}
    isempty(attractor) && throw(ArgumentError("attractor collection must not be empty"))
    sources = CartesianIndex{D}[]
    sizehint!(sources, length(attractor))
    for (i, state) in pairs(attractor)
        length(state) == D || throw(
            DimensionMismatch(
                "attractor state $i has length $(length(state)) but sys has D=$D",
            ),
        )
        all(x -> x isa Real, state) || throw(
            ArgumentError("attractor state $i must contain only real coordinates"),
        )
        push!(sources, _source_cell(grid, SVector{D, T}(state)))
    end
    unique!(sources)
    return sources
end

function _olim_periodic_axes(bc, ::Val{D}) where {D}
    bcs = _normalize_bc(bc, Val(D))
    @inbounds for k in 1:D
        bcs[k] isa Absorbing && throw(
            ArgumentError(
                "quasipotential supports Reflecting and Periodic boundaries; " *
                    "Absorbing is a generator boundary condition and does not define " *
                    "an OLIM boundary value",
            ),
        )
    end
    return ntuple(k -> bcs[k] isa Periodic, Val(D))
end

"""
    quasipotential(sys, grid, attractor; band_radius, near_source_layers,
                   bc, verbose, show_progress) -> QuasiPotential{D}

Compute the Freidlin-Wentzell quasipotential field `U_A(x)` from `attractor`
using the Ordered Line Integral Method (Dahiya and Cameron 2018). The state
dimension `D` is taken from `sys::CoupledSDEs{IIP, D, I, P}` and must match
`grid::CartesianGrid{D}`. A warning is emitted for `D > 4`.

`attractor` may be one state vector (a fixed point), or a vector of state
vectors describing an extended attractor such as a limit cycle. Every grid
cell touched by an extended attractor is initialised with the known value
`U = 0`; no derivative or curvature data are used. Analytic CARE seeding is
only available for a one-state attractor, so `near_source_layers` defaults to
`0` for an extended attractor and must remain zero.

# Keyword arguments
* `band_radius::Int  = default_K(grid)`: accepted-band radius in grid cells.
* `near_source_layers::Int = 3`: size of the analytic CARE seed box;
  `0` disables analytic seeding.
* `bc = Reflecting()`: grid topology, either one boundary condition for every
  axis or a `D`-tuple for per-axis control. [`Reflecting`](@ref) keeps the OLIM
  stencil inside the grid; [`Periodic`](@ref) wraps neighbours and line-integral
  geometry across the corresponding seam. [`Absorbing`](@ref) is rejected
  because it does not specify a Hamilton-Jacobi boundary value.
* `regularization::Real = default_regularization(grid)`: for a system with a single
  noiseless coordinate (rank-1 diffusion, e.g. momentum-only noise in a second-order
  Langevin or van der Pol oscillator) the trace-normalized diffusion is singular and
  has no metric. This amount is added to that coordinate's diagonal (a structural
  vanishing-viscosity perturbation that leaves the noisy block exact) to make it
  invertible; the degenerate Freidlin-Wentzell quasipotential is recovered as it tends
  to 0. The default scales like `1 / minimum(grid.nbox)` (about `0.04` at 80 cells), so
  it shrinks under refinement; smaller values are more accurate but stiffer. The escape
  sheet and saddle barrier carry negligible bias, the non-escape region carries the
  regularization bias. Ignored for full-rank systems; two or more noiseless coordinates
  (sub-Riemannian) or a non-coordinate-aligned null space are rejected.
* `verbose::Bool      = false`
* `show_progress::Bool = true`
* `threaded::Bool = true`: spread each pop's local updates over the available Julia
  threads. The result is *bitwise* identical to the serial sweep (see `_sweep_loop!`),
  so this trades nothing for speed and is on by default; it is a no-op when Julia runs
  with one thread. Pass `false` when `quasipotential` is itself being called from
  inside a parallel region and you would rather own the threads yourself.
* `_full_rescan::Bool = false`: internal. Rebuilds each candidate's value from every
  front simplex in its stencil on every pop, instead of only from the simplexes the
  newly popped cell just created (see `_candidate_update`). Measured 27x slower at
  `41x41` with `K = 8` and 64x at `61x61` with `K = 12` -- the gap widens with the
  stencil, since that is exactly the work being skipped. Present only so the tests can
  pin the two against each other.

# Accuracy

The sweep minimises over a restricted set of paths (polylines through the accepted
front), so its error is one-sided: `U` is an *upper* estimate, and `U >= 0` always.
A cell below the exact value is a bug, not discretisation.

How large that overshoot is depends almost entirely on how anisotropic the metric is:

* Full-rank, well-conditioned diffusion: essentially exact. For `b = -x` with `a = I`
  the error is `0.13%` at worst at `31` cells across and `0.03%` at `121`, with a
  median near `10^-5`.
* Rank-1 diffusion with `regularization = δ`: the metric has anisotropy
  `sqrt(2/δ)`, and the error grows as `δ` shrinks. On a linear oscillator against an
  exact CARE reference, over the inner 60% of the box (outside that the true optimal
  path leaves the domain, which is truncation rather than method error):

  | cells | `δ = 0.05` median / worst | `δ = 0.005` median / worst |
  |-------|---------------------------|----------------------------|
  |  61   |      5.9% / 32%           |       9.2% / 46%           |
  | 121   |      3.6% / 18%           |       5.8% / 23%           |

  i.e. falling like `O(h)`, and positive in every cell, as an upper bound must be.
* `band_radius` past the point where the stencil already contains the cheapest path
  changes the interior not at all: at `61` cells and `δ = 0.005`, `K = 8` through `24`
  agree to the last bit over the inner 60%, and `K = 5` differs there by only `7e-3`.
  What a wider stencil does buy is the *outer* region, where it keeps improving well
  past `default_K`. So lower `band_radius` freely if you only care about the interior,
  and raise it if you care about `U` near the domain edge.
* Second derivatives of `U` are far less accurate than `U`. Differencing a field
  known to `O(h)` twice leaves errors of tens of percent on a transverse Hessian even
  where `U` itself is good to a fraction of a percent; fit `U` over a finite window
  and account for its non-quadratic terms rather than differencing pointwise.

# Performance

The sweep does `O(N)` heap pops, each updating the `O(K^D)` unfinalised cells in the
popped cell's stencil. Each such update costs `O(1)`, not `O(K^D)`: it only examines the
simplexes the newly popped cell just created, because every other front simplex was
already accounted for at an earlier pop (see `_candidate_update`). The work is dominated
by evaluating `Φ`, a three-node Simpson quadrature of the Lagrangian, a few times per
candidate simplex, so an expensive `dynamic_rule` is felt directly.

In practice, for a 2D rank-1 problem:

| grid  | `K`  | serial  | threaded (6 cores) |
|-------|------|---------|--------------------|
| 41x41 |   8  | 0.047 s |      0.035 s       |
| 61x61 |  12  | 0.22 s  |      0.12 s        |
| 81x81 |   9  | 0.25 s  |      0.14 s        |

* Refining the grid is the main cost: `O(N)` pops times `O(K^D)` candidates, so 2x in
  each 2D direction is about 4x, plus whatever `default_K` adds by growing like
  `sqrt(n)`.
* `band_radius` is much weaker than it looks. Measured cost grows like `K^1.6` at fixed
  grid (`K = 5 -> 24` is 11x), not like the stencil volume, because most cells in a wide
  stencil are already finalised and skipped.
* Threading is on by default and exact; see `threaded` above. Its benefit is modest
  (`1.3-1.9x` here) because what remains after the parallel batch -- the heap and the
  serial apply pass -- is now a real fraction of the total.

See also: [`QuasiPotential`](@ref), [`BackRef`](@ref).
"""
function quasipotential(
        sys::CoupledSDEs{IIP, D, I, P},
        grid::CartesianGrid{D, T},
        attractor::AbstractVector{<:Real};
        band_radius::Int = default_K(grid),
        near_source_layers::Int = 3,
        regularization::Real = default_regularization(grid),
        bc = Reflecting(),
        verbose::Bool = false,
        show_progress::Bool = true,
        threaded::Bool = true,
        _full_rescan::Bool = false,
    ) where {IIP, D, I, P, T}
    length(attractor) == D || throw(
        DimensionMismatch(
            "attractor has length $(length(attractor)) but sys has D=$D",
        ),
    )
    proper_FW_system(sys)
    D > 4 && @warn "quasipotential in D=$D: per-axis grid resolution will be coarse" maxlog = 1
    periodic = _olim_periodic_axes(bc, Val(D))
    x_A = SVector{D, T}(attractor)
    return _quasipotential_impl(
        sys, grid, x_A,
        Val(band_radius), Val(near_source_layers),
        T(regularization), periodic, verbose, show_progress, threaded, _full_rescan,
    )
end

function quasipotential(
        sys::CoupledSDEs{IIP, D, I, P},
        grid::CartesianGrid{D, T},
        attractor::AbstractVector{<:AbstractVector};
        band_radius::Int = default_K(grid),
        near_source_layers::Int = 0,
        regularization::Real = default_regularization(grid),
        bc = Reflecting(),
        verbose::Bool = false,
        show_progress::Bool = true,
        threaded::Bool = true,
        _full_rescan::Bool = false,
    ) where {IIP, D, I, P, T}
    near_source_layers == 0 || throw(
        ArgumentError(
            "near_source_layers must be 0 for an extended attractor; " *
                "analytic CARE seeding is only defined at a stable fixed point",
        ),
    )
    proper_FW_system(sys)
    D > 4 && @warn "quasipotential in D=$D: per-axis grid resolution will be coarse" maxlog = 1
    sources = _source_cells(grid, attractor)
    periodic = _olim_periodic_axes(bc, Val(D))
    return _quasipotential_impl(
        sys, grid, sources,
        Val(band_radius), T(regularization), periodic,
        verbose, show_progress, threaded, _full_rescan,
    )
end

function _quasipotential_impl(
        sys::CoupledSDEs{IIP, D, I, P},
        grid::CartesianGrid{D, T},
        x_A::SVector{D, T},
        ::Val{K}, ::Val{K_seed},
        regularization::T, periodic::NTuple{D, Bool},
        verbose::Bool, show_progress::Bool, threaded::Bool, full_rescan::Bool,
    ) where {IIP, D, I, P, T, K, K_seed}
    src = _source_cell(grid, x_A)
    L = _geometric_lagrangian(sys, T; regularization = regularization)
    state = _OLIMState(grid, T, L.Q_inv isa SMatrix; periodic = periodic)
    _sweep!(
        state, grid, src, sys, L, Val(K), Val(K_seed);
        verbose = verbose, show_progress = show_progress, threaded = threaded,
        full_rescan = full_rescan,
    )
    return QuasiPotential{D, T}(state.U, state.back_pointer, src, grid)
end

function _quasipotential_impl(
        sys::CoupledSDEs{IIP, D, I, P},
        grid::CartesianGrid{D, T},
        sources::Vector{CartesianIndex{D}},
        ::Val{K}, regularization::T, periodic::NTuple{D, Bool},
        verbose::Bool, show_progress::Bool, threaded::Bool, full_rescan::Bool,
    ) where {IIP, D, I, P, T, K}
    L = _geometric_lagrangian(sys, T; regularization = regularization)
    state = _OLIMState(grid, T, L.Q_inv isa SMatrix; periodic = periodic)
    @inbounds for source in sources
        state.U[source] = zero(T)
        state.status[source] = _ACCEPTED
        state.back_pointer[source] = BackRef{D}(source, source, NaN32)
    end
    src = first(sources)
    _sweep!(
        state, grid, src, sys, L, Val(K), Val(0);
        verbose = verbose, show_progress = show_progress, threaded = threaded,
        full_rescan = full_rescan,
    )
    return QuasiPotential{D, T}(
        state.U, state.back_pointer, src, grid, sources,
    )
end
