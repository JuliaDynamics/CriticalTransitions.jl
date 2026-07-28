# Per-cell drift cache (additive noise only): one drift eval per cell up front
# replaces the thousands of redundant cell-center evals across the sweep.
@inline function _fill_node_cache!(
        state::_OLIMState{D, T}, grid::CartesianGrid{D, T},
        L::_GeometricLagrangian{D, T, B, A},
    ) where {D, T, B, A <: SMatrix}
    Qinv = L.Q_inv
    @inbounds for I in CartesianIndices(state.U)
        b = L.b(cell_center(grid, I))
        q = dot(b, Qinv * b)
        state.b_c[I] = b
        state.q_c[I] = q
        state.sq_c[I] = sqrt(q)
    end
    return state
end
@inline _fill_node_cache!(
    state::_OLIMState{D, T}, ::CartesianGrid{D, T}, ::_GeometricLagrangian{D, T},
) where {D, T} = state

# Cached drift data at cell `I` (additive noise) or `nothing` (multiplicative).
@inline function _node_at(
        state::_OLIMState{D, T}, I::CartesianIndex{D},
        ::_GeometricLagrangian{D, T, B, A},
    ) where {D, T, B, A <: SMatrix}
    @inbounds return _NodeData{D, T}(state.b_c[I], state.q_c[I], state.sq_c[I])
end
@inline _node_at(
    ::_OLIMState{D, T}, ::CartesianIndex{D}, ::_GeometricLagrangian{D, T},
) where {D, T} = nothing

@inline function _vertex_candidate(
        c_x::SVector{D, T}, c_y::SVector{D, T},
        U_y::T, L::_GeometricLagrangian{D, T}, nd_y, nd_x,
    ) where {D, T}
    v = c_x - c_y
    return U_y + _line_integral(L, c_y, v, nd_y, nd_x)
end

# The interior point `y` is interpolated (never a cell center), so the s=0 node is
# always live; only the shared s=1 node `c_x` is cached via `nd_x`.
#
# `s0`/`s1` are *prepared* slopes (secant-substituted and limited by
# `_prepare_slopes`). They do not depend on `λ`, so the caller hoists that work out
# of the many evaluations made along a single edge.
@inline function _edge_phi(
        c_x::SVector{D, T}, c_y0::SVector{D, T}, c_y1::SVector{D, T},
        U0::T, U1::T, s0::T, s1::T, λ::T,
        L::_GeometricLagrangian{D, T}, nd_x,
    ) where {D, T}
    y = (one(T) - λ) * c_y0 + λ * c_y1
    return _hermite_prepared(U0, U1, s0, s1, λ) +
        _line_integral(L, y, c_x - y, nothing, nd_x)
end

@inline function _edge_dphi(
        c_x::SVector{D, T}, c_y0::SVector{D, T}, c_y1::SVector{D, T},
        U0::T, U1::T, s0::T, s1::T, λ::T,
        L::_GeometricLagrangian{D, T}, nd_x,
    ) where {D, T}
    h = sqrt(eps(T))
    λp = clamp(λ + h, zero(T), one(T))
    λm = clamp(λ - h, zero(T), one(T))
    return (
        _edge_phi(c_x, c_y0, c_y1, U0, U1, s0, s1, λp, L, nd_x) -
            _edge_phi(c_x, c_y0, c_y1, U0, U1, s0, s1, λm, L, nd_x)
    ) / (λp - λm)
end

"""
Convergence tolerance on the edge parameter `λ` for the interior root find, and the
main lever on how often `Φ` is evaluated: `itp_root` takes on the order of
`log2(1 / (2 * _EDGE_ATOL))` iterations at two `Φ` evaluations each.

It can afford to be loose, because `Φ` is stationary at `λ*`: mislocating the minimiser
by `ε` perturbs `U` only by `Φ'' ε^2 / 2`, and erring high is the safe direction for an
upper-bound method. Loosening it does *not* risk missing the narrow interior valley an
anisotropic metric produces -- a bracketed root find converges to the root of `dΦ/dλ`
whatever the tolerance, and only its precision changes. Measured on a rank-1 oscillator
at `61x61`, going from `1e-4` to `1e-2` moved the median error in `U` by 0.02 percentage
points, so this is not where the accuracy lives; it is kept tight because since
`_candidate_update` went incremental the root find is no longer the bottleneck either.
"""
const _EDGE_ATOL = 1.0e-4

"""
Edge candidate minimum (§4.4). Evaluates Φ(0), Φ(1), and an ITP root of
dΦ/dλ; the smallest wins. `λ* ∈ {0, 1}` indicates a boundary win.
"""
struct _EdgeDPhi{D, T, LT, NX}
    c_x::SVector{D, T}
    c_y0::SVector{D, T}
    c_y1::SVector{D, T}
    U0::T
    U1::T
    s0::T
    s1::T
    L::LT
    nd_x::NX
end
@inline (g::_EdgeDPhi)(λ) =
    _edge_dphi(g.c_x, g.c_y0, g.c_y1, g.U0, g.U1, g.s0, g.s1, λ, g.L, g.nd_x)

# `Φ0`/`Φ1` are the edge's λ=0 and λ=1 values. They coincide with the *vertex*
# candidates from `c_y0` and `c_y1`, which every caller has already computed, so they
# are required rather than defaulted -- recomputing them here would double the cheap
# part of the edge and hide it behind a default argument.
@inline function _edge_minimum(
        c_x::SVector{D, T},
        c_y0::SVector{D, T}, c_y1::SVector{D, T},
        U0::T, U1::T, m0::T, m1::T,
        L::_GeometricLagrangian{D, T}, nd_x, Φ0::T, Φ1::T,
    ) where {D, T}
    s0, s1 = _prepare_slopes(U0, U1, m0, m1)
    best, λstar = (Φ0 <= Φ1) ? (Φ0, zero(T)) : (Φ1, one(T))
    # The interior root is always sought. It used to be skipped when both endpoints
    # exceeded the caller's running best, on the assumption that the interior minimum
    # is bracketed by smaller endpoints; it is not. Under an anisotropic metric Φ(λ)
    # has a narrow interior valley (the foot point where the segment aligns with the
    # drift), so the endpoints are a poor proxy and the only sound bound available,
    # Φ ≥ min(U0, U1), is far too weak to prune on. Skipping cost up to a factor 1.8
    # in worst-case accuracy and made the answer depend on cell visit order.
    #
    # Entering the root find used to cost four Φ evaluations: `itp_root` opens with
    # `f(0)` and `f(1)`, and `f = _edge_dphi` is a two-point difference. But
    # `_edge_dphi` clamps its stencil at the boundary, so its endpoint values are
    # *one-sided* differences whose far node is Φ(0) resp. Φ(1) -- both already in
    # hand. Building them here costs two Φ evaluations instead of four and is
    # bitwise identical (the clamped denominators are exactly `h`).
    h = sqrt(eps(T))
    λ1m = one(T) - h
    dΦ0 = (_edge_phi(c_x, c_y0, c_y1, U0, U1, s0, s1, h, L, nd_x) - Φ0) / h
    dΦ1 = (Φ1 - _edge_phi(c_x, c_y0, c_y1, U0, U1, s0, s1, λ1m, L, nd_x)) / (one(T) - λ1m)
    # Only a bracketed sign change can yield an interior minimiser. Without one
    # `itp_root` either reports `ok = false` or returns an endpoint, and an endpoint
    # `λroot` fails the `0 < λroot < 1` test below -- so this gate skips exactly the
    # calls whose result is discarded, and skips them before the setup cost.
    dΦ0 * dΦ1 < zero(T) || return (best, λstar)
    g = _EdgeDPhi(c_x, c_y0, c_y1, U0, U1, s0, s1, L, nd_x)
    λroot, ok = itp_root(
        g, zero(T), one(T); atol = T(_EDGE_ATOL), maxiter = 20, fa = dΦ0, fb = dΦ1,
    )
    if ok && zero(T) < λroot < one(T)
        Φi = _edge_phi(c_x, c_y0, c_y1, U0, U1, s0, s1, λroot, L, nd_x)
        if Φi < best
            best, λstar = Φi, λroot
        end
    end
    return (best, λstar)
end

# A cell's `U` is final once it has been popped from the heap (`_FRONT`) or pruned
# (`_ACCEPTED`). A `_CONSIDERED` cell still holds a tentative upper bound that later
# updates can only lower, so a divided difference taken against it is biased high --
# and a slope that is too large is exactly what drives the Hermite interpolant below
# both edge endpoints. `_UNKNOWN` cells hold `Inf` and are excluded by `isfinite`.
@inline function _is_final(state::_OLIMState{D}, I::CartesianIndex{D}) where {D}
    @inbounds s = state.status[I]
    return s == _FRONT || s == _ACCEPTED
end

@inline function _edge_slope_estimates(
        state::_OLIMState{D, T},
        y0::CartesianIndex{D}, y1::CartesianIndex{D},
    ) where {D, T}
    dir = y1 - y0
    nbox = state.nbox
    y0m = _shift_in_bounds(y0, -dir, nbox)
    y1p = _shift_in_bounds(y1, dir, nbox)
    usable(I) = I !== nothing && _is_final(state, I) && isfinite(state.U[I])
    m0 = usable(y0m) ? (state.U[y1] - state.U[y0m]) / 2 : T(NaN)
    m1 = usable(y1p) ? (state.U[y1p] - state.U[y0]) / 2 : T(NaN)
    return (m0, m1)
end

"""
Active-set minimisation on Δ_2: evaluate 3 vertices and 3 edges (ITP root),
then attempt an interior Newton step with FD gradient/Hessian.
"""
@inline function _tri_phi(
        c_x::SVector{3, T},
        c_y0::SVector{3, T}, c_y1::SVector{3, T}, c_y2::SVector{3, T},
        U0::T, U1::T, U2::T,
        L::_GeometricLagrangian{3, T}, λ1::T, λ2::T, nd_x,
    ) where {T}
    y = (one(T) - λ1 - λ2) * c_y0 + λ1 * c_y1 + λ2 * c_y2
    return ((one(T) - λ1 - λ2) * U0 + λ1 * U1 + λ2 * U2) +
        _line_integral(L, y, c_x - y, nothing, nd_x)
end

@inline function _triangle_minimum(
        c_x::SVector{3, T},
        c_y0::SVector{3, T}, c_y1::SVector{3, T}, c_y2::SVector{3, T},
        U0::T, U1::T, U2::T,
        L::_GeometricLagrangian{3, T}, nd_x,
    ) where {T}
    Φ00 = _tri_phi(c_x, c_y0, c_y1, c_y2, U0, U1, U2, L, zero(T), zero(T), nd_x)
    best, bλ1, bλ2 = Φ00, zero(T), zero(T)
    Φ1 = _tri_phi(c_x, c_y0, c_y1, c_y2, U0, U1, U2, L, one(T), zero(T), nd_x)
    Φ1 < best && ((best, bλ1, bλ2) = (Φ1, one(T), zero(T)))
    Φ2 = _tri_phi(c_x, c_y0, c_y1, c_y2, U0, U1, U2, L, zero(T), one(T), nd_x)
    Φ2 < best && ((best, bλ1, bλ2) = (Φ2, zero(T), one(T)))

    # The three barycentric corners are exactly the edges' λ=0 / λ=1 values, so
    # Φ00/Φ1/Φ2 above serve as the endpoint values for all three edges.
    Φe01, λs01 = _edge_minimum(
        c_x, c_y0, c_y1, U0, U1, T(NaN), T(NaN), L, nd_x, Φ00, Φ1,
    )
    if Φe01 < best
        best = Φe01; bλ1 = λs01; bλ2 = zero(T)
    end
    Φe02, λs02 = _edge_minimum(
        c_x, c_y0, c_y2, U0, U2, T(NaN), T(NaN), L, nd_x, Φ00, Φ2,
    )
    if Φe02 < best
        best = Φe02; bλ1 = zero(T); bλ2 = λs02
    end
    Φe12, λs12 = _edge_minimum(
        c_x, c_y1, c_y2, U1, U2, T(NaN), T(NaN), L, nd_x, Φ1, Φ2,
    )
    if Φe12 < best
        best = Φe12; bλ1 = one(T) - λs12; bλ2 = λs12
    end

    h = sqrt(eps(T))
    λ1 = one(T) / 3; λ2 = one(T) / 3
    @inline P(a, b) = _tri_phi(c_x, c_y0, c_y1, c_y2, U0, U1, U2, L, a, b, nd_x)
    for _ in 1:8
        gx = (P(λ1 + h, λ2) - P(λ1 - h, λ2)) / (2h)
        gy = (P(λ1, λ2 + h) - P(λ1, λ2 - h)) / (2h)
        max(abs(gx), abs(gy)) < T(1.0e-9) && break
        Hxx = (P(λ1 + h, λ2) - 2P(λ1, λ2) + P(λ1 - h, λ2)) / h^2
        Hyy = (P(λ1, λ2 + h) - 2P(λ1, λ2) + P(λ1, λ2 - h)) / h^2
        Hxy = (
            P(λ1 + h, λ2 + h) - P(λ1 + h, λ2 - h) -
                P(λ1 - h, λ2 + h) + P(λ1 - h, λ2 - h)
        ) / (4h^2)
        det = Hxx * Hyy - Hxy^2
        abs(det) < eps(T) && break
        dλ1 = (Hyy * gx - Hxy * gy) / det
        dλ2 = (-Hxy * gx + Hxx * gy) / det
        λ1 -= dλ1; λ2 -= dλ2
        (λ1 < 0 || λ2 < 0 || λ1 + λ2 > 1) && break
    end
    if λ1 > 0 && λ2 > 0 && λ1 + λ2 < 1
        Φi = _tri_phi(c_x, c_y0, c_y1, c_y2, U0, U1, U2, L, λ1, λ2, nd_x)
        Φi < best && ((best, bλ1, bλ2) = (Φi, λ1, λ2))
    end
    return (best, bλ1, bλ2)
end

"""
Iterates the K-stencil and considers vertex + simplex (D=2 edge, D=3
triangle) candidates. Returns `(Inf, BackRef{D}())` if no usable neighbour.
"""
function _local_update(
        ::Val{D}, x::CartesianIndex{D},
        state::_OLIMState{D, T},
        grid::CartesianGrid{D, T},
        L::_GeometricLagrangian{D, T},
        ::Val{K},
    ) where {D, T, K}
    nbox = state.nbox
    c_x = cell_center(grid, x)
    nd_x = _node_at(state, x, L)
    best = T(Inf)
    best_ref = BackRef{D}()
    offsets = _stencil_offsets(Val(K), Val(D))

    # Vertex candidates from ACCEPTED-only cells (FRONT cells are handled
    # together with their edges in the simplex pass to share Φ(0) work).
    @inbounds for δ in offsets
        y = _shift_in_bounds(x, δ, nbox); y === nothing && continue
        state.status[y] == _ACCEPTED || continue
        U_y = state.U[y]
        (isfinite(U_y) && U_y < best) || continue
        c_y = cell_center(grid, y)
        Φ = _vertex_candidate(c_x, c_y, U_y, L, _node_at(state, y, L), nd_x)
        if Φ < best
            best = Φ
            best_ref = BackRef{D}(y, y, NaN32)
        end
    end
    best, best_ref = _add_simplex_candidates(
        best, best_ref, x, c_x, state, grid, L, Val(K), Val(D),
    )
    return (best, best_ref)
end

@inline _add_simplex_candidates(
    best, best_ref, x, c_x, state, grid, L, ::Val{K}, ::Val{D},
) where {K, D} = (best, best_ref)

"""
Update `x` in response to a single cell `c` having just been popped.

Popping `c` is the only status change, so the only candidates `x` has not already been
offered are the ones built from `c`: the vertex `c` itself and the simplexes having `c`
as a corner. Every other front simplex was already complete at an earlier pop -- the
stencil relation `|x - y| <= K` is symmetric, so when a simplex's last corner `y` was
popped, `x` lay in `y`'s stencil exactly when `y` lies in `x`'s, and `x` was updated
then. `state.U[x]` therefore already holds the minimum over all of them, and this starts
from that value instead of rebuilding it.

That turns the per-pop cost for a cell from "every front simplex in the `K`-disc"
(hundreds of edges at the radii an anisotropic metric needs) into "the at most eight
edges touching `c`", which is where nearly all of the sweep's time was going.

The one thing the full rescan did additionally was re-evaluate old edges whose slope
estimates had since improved, as more neighbours became final. Skipping that can only
leave `U` higher, never lower, so the sweep stays the upper bound it is documented to be.
The price is small: on a rank-1 oscillator at `61x61` with `K = 12`, against an exact CARE
reference, the median error goes from 9.07% to 9.20% and the worst case is unchanged at
46%, with `U` raised rather than lowered in 99.9% of cells. `_full_rescan = true` on
`quasipotential` restores the old behaviour, and the tests hold the two side by side.
"""
@inline _candidate_update(
    v::Val{D}, x, _c, state, grid, L, vk, _full,
) where {D} = _local_update(v, x, state, grid, L, vk)

# Only `D = 2` and `D = 3` may take the incremental path. For every other `D` the
# simplex pass is a no-op, so candidates come solely from the `_ACCEPTED` vertex loop,
# and a cell joins that set by being *pruned* (`_FRONT -> _ACCEPTED`) -- an event that
# triggers no update and so cannot be caught incrementally.
for VD in (:(Val{2}), :(Val{3}))
    @eval @inline function _candidate_update(v::$VD, x, c, state, grid, L, vk, full)
        # A cell seen for the first time has no accumulated value to build on.
        @inbounds (full || state.status[x] == _UNKNOWN) &&
            return _local_update(v, x, state, grid, L, vk)
        return _incremental_update(v, x, c, state, grid, L, vk)
    end
end

function _incremental_update(
        ::Val{2}, x::CartesianIndex{2}, c::CartesianIndex{2},
        state::_OLIMState{2, T}, grid::CartesianGrid{2, T},
        L::_GeometricLagrangian{2, T}, ::Val{K},
    ) where {T, K}
    nbox = state.nbox
    front = state.front
    @inbounds best = state.U[x]
    @inbounds best_ref = state.back_pointer[x]
    @inbounds U_c = state.U[c]
    isfinite(U_c) || return (best, best_ref)
    c_x = cell_center(grid, x)
    nd_x = _node_at(state, x, L)
    c_c = cell_center(grid, c)
    nd_c = _node_at(state, c, L)
    Φc = T(NaN)
    if U_c < best
        Φc = U_c + _line_integral(L, c_c, c_x - c_c, nd_c, nd_x)
        if Φc < best
            best = Φc
            best_ref = BackRef{2}(c, c, NaN32)
        end
    end
    @inbounds for δ in _chebyshev_neighbors(Val(2))
        y = _shift_in_bounds(c, δ, nbox)
        y === nothing && continue
        front[y] || continue
        _in_circular_stencil(x, y, Val(K)) || continue
        U_y = state.U[y]
        isfinite(U_y) || continue
        min(U_c, U_y) < best || continue
        isnan(Φc) && (Φc = U_c + _line_integral(L, c_c, c_x - c_c, nd_c, nd_x))
        c_y = cell_center(grid, y)
        Φy = U_y + _line_integral(L, c_y, c_x - c_y, _node_at(state, y, L), nd_x)
        # The full scan always orients an edge from its lower linear index. Match that,
        # so the interpolant and the slope stencil are the same whichever end `c` is.
        lo = c < y
        y0, y1 = lo ? (c, y) : (y, c)
        cy0, cy1 = lo ? (c_c, c_y) : (c_y, c_c)
        U0, U1 = lo ? (U_c, U_y) : (U_y, U_c)
        Φ0, Φ1 = lo ? (Φc, Φy) : (Φy, Φc)
        m0, m1 = _edge_slope_estimates(state, y0, y1)
        Φ, λstar = _edge_minimum(c_x, cy0, cy1, U0, U1, m0, m1, L, nd_x, Φ0, Φ1)
        if Φ < best
            best = Φ
            best_ref = BackRef{2}(y0, y1, Float32(λstar))
        end
    end
    return (best, best_ref)
end

function _incremental_update(
        ::Val{3}, x::CartesianIndex{3}, c::CartesianIndex{3},
        state::_OLIMState{3, T}, grid::CartesianGrid{3, T},
        L::_GeometricLagrangian{3, T}, ::Val{K},
    ) where {T, K}
    nbox = state.nbox
    front = state.front
    @inbounds best = state.U[x]
    @inbounds best_ref = state.back_pointer[x]
    @inbounds U_c = state.U[c]
    isfinite(U_c) || return (best, best_ref)
    c_x = cell_center(grid, x)
    nd_x = _node_at(state, x, L)
    c_c = cell_center(grid, c)
    Φc = U_c + _line_integral(L, c_c, c_x - c_c, _node_at(state, c, L), nd_x)
    if Φc < best
        best = Φc
        best_ref = BackRef{3}(c, c, NaN32)
    end
    cheb = _chebyshev_neighbors(Val(3))
    # Every corner of a triangle the full scan enumerates is Chebyshev-adjacent to the
    # other two, so the triangles containing `c` are exactly the mutually adjacent
    # front pairs drawn from `c`'s own neighbourhood.
    @inbounds for δ1 in cheb
        p = _shift_in_bounds(c, δ1, nbox)
        p === nothing && continue
        front[p] || continue
        _in_circular_stencil(x, p, Val(K)) || continue
        isfinite(state.U[p]) || continue
        for δ2 in cheb
            q = _shift_in_bounds(c, δ2, nbox)
            q === nothing && continue
            q > p || continue
            front[q] || continue
            _is_chebyshev_adjacent(p, q) || continue
            _in_circular_stencil(x, q, Val(K)) || continue
            isfinite(state.U[q]) || continue
            # Sort the corners so the barycentric frame matches the full scan's, which
            # always bases a triangle at its lowest-index corner.
            a, b, d = _sorted3(c, p, q)
            Φ, λ1, _ = _triangle_minimum(
                c_x, cell_center(grid, a), cell_center(grid, b), cell_center(grid, d),
                state.U[a], state.U[b], state.U[d], L, nd_x,
            )
            if Φ < best
                best = Φ
                best_ref = BackRef{3}(a, b, Float32(λ1))
            end
        end
    end
    return (best, best_ref)
end

@inline function _sorted3(p::T, q::T, r::T) where {T}
    lo, hi = q < r ? (q, r) : (r, q)
    p < lo && return (p, lo, hi)
    p < hi && return (lo, p, hi)
    return (lo, hi, p)
end

function _add_simplex_candidates(
        best, best_ref, x::CartesianIndex{2},
        c_x::SVector{2, T},
        state::_OLIMState{2, T},
        grid::CartesianGrid{2, T},
        L::_GeometricLagrangian{2, T},
        ::Val{K}, ::Val{2},
    ) where {T, K}
    nbox = state.nbox
    front = state.front
    nd_x = _node_at(state, x, L)
    offsets = _stencil_offsets(Val(K), Val(2))
    cheb = _chebyshev_neighbors(Val(2))
    @inbounds for δ0 in offsets
        y0 = _shift_in_bounds(x, δ0, nbox)
        y0 === nothing && continue
        front[y0] || continue
        U0 = state.U[y0]
        isfinite(U0) || continue
        c_y0 = cell_center(grid, y0)
        nd_y0 = _node_at(state, y0, L)
        # Φ0 = vertex candidate from y0 (same as edge Φ(λ=0)). Computed at most once
        # per y0, then reused for the vertex update and as the Φ(0) endpoint of every
        # edge leaving y0. Deferred because `U0 ≥ best` bounds Φ0 out on its own.
        Φ0 = T(NaN)
        if U0 < best
            Φ0 = U0 + _line_integral(L, c_y0, c_x - c_y0, nd_y0, nd_x)
            if Φ0 < best
                best = Φ0
                best_ref = BackRef{2}(y0, y0, NaN32)
            end
        end
        for δch in cheb
            y1 = _shift_in_bounds(y0, δch, nbox)
            y1 === nothing && continue
            y1 > y0 || continue
            front[y1] || continue
            _in_circular_stencil(x, y1, Val(K)) || continue
            U1 = state.U[y1]
            isfinite(U1) || continue
            # The limited interpolant keeps the edge value in [min(U0, U1), max(U0, U1)]
            # and the line integral is nonnegative, so Φ ≥ min(U0, U1). Requiring *both*
            # endpoints below `best` (as this once did) discards edges whose cheaper end
            # can still win, and since each edge is visited only from its lower-index
            # endpoint those candidates were lost outright.
            min(U0, U1) < best || continue
            isnan(Φ0) && (Φ0 = U0 + _line_integral(L, c_y0, c_x - c_y0, nd_y0, nd_x))
            c_y1 = cell_center(grid, y1)
            # Φ(1) is the vertex candidate from y1, reused as the edge's λ=1 endpoint.
            Φ1 = U1 + _line_integral(L, c_y1, c_x - c_y1, _node_at(state, y1, L), nd_x)
            m0, m1 = _edge_slope_estimates(state, y0, y1)
            Φ, λstar = _edge_minimum(
                c_x, c_y0, c_y1, U0, U1, m0, m1, L, nd_x, Φ0, Φ1,
            )
            if Φ < best
                best = Φ
                best_ref = BackRef{2}(y0, y1, Float32(λstar))
            end
        end
    end
    return (best, best_ref)
end

function _add_simplex_candidates(
        best, best_ref, x::CartesianIndex{3},
        c_x::SVector{3, T},
        state::_OLIMState{3, T},
        grid::CartesianGrid{3, T},
        L::_GeometricLagrangian{3, T},
        ::Val{K}, ::Val{3},
    ) where {T, K}
    nbox = state.nbox
    front = state.front
    nd_x = _node_at(state, x, L)
    offsets = _stencil_offsets(Val(K), Val(3))
    cheb = _chebyshev_neighbors(Val(3))
    @inbounds for δ0 in offsets
        y0 = _shift_in_bounds(x, δ0, nbox)
        y0 === nothing && continue
        front[y0] || continue
        U0 = state.U[y0]
        isfinite(U0) || continue
        c_y0 = cell_center(grid, y0)
        Φ0 = U0 + _line_integral(L, c_y0, c_x - c_y0, _node_at(state, y0, L), nd_x)
        if Φ0 < best
            best = Φ0
            best_ref = BackRef{3}(y0, y0, NaN32)
        end
        for δch1 in cheb
            y1 = _shift_in_bounds(y0, δch1, nbox)
            y1 === nothing && continue
            y1 > y0 || continue
            front[y1] || continue
            _in_circular_stencil(x, y1, Val(K)) || continue
            U1 = state.U[y1]
            isfinite(U1) || continue
            c_y1 = cell_center(grid, y1)
            for δch2 in cheb
                y2 = _shift_in_bounds(y0, δch2, nbox)
                y2 === nothing && continue
                y2 > y1 || continue
                front[y2] || continue
                _is_chebyshev_adjacent(y1, y2) || continue
                _in_circular_stencil(x, y2, Val(K)) || continue
                U2 = state.U[y2]
                isfinite(U2) || continue
                c_y2 = cell_center(grid, y2)
                Φ, λ1, _ = _triangle_minimum(c_x, c_y0, c_y1, c_y2, U0, U1, U2, L, nd_x)
                if Φ < best
                    best = Φ
                    best_ref = BackRef{3}(y0, y1, Float32(λ1))
                end
            end
        end
    end
    return (best, best_ref)
end
