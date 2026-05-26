@inline function _vertex_candidate(
        c_x::SVector{D, T}, c_y::SVector{D, T},
        U_y::T, L::_GeometricLagrangian{D, T},
    ) where {D, T}
    v = c_x - c_y
    return U_y + _line_integral(L, c_y, v)
end

@inline function _edge_phi(
        c_x::SVector{D, T}, c_y0::SVector{D, T}, c_y1::SVector{D, T},
        U0::T, U1::T, m0::T, m1::T, λ::T,
        L::_GeometricLagrangian{D, T},
    ) where {D, T}
    y = (one(T) - λ) * c_y0 + λ * c_y1
    return _hermite_U(U0, U1, m0, m1, λ) + _line_integral(L, y, c_x - y)
end

@inline function _edge_dphi(
        c_x::SVector{D, T}, c_y0::SVector{D, T}, c_y1::SVector{D, T},
        U0::T, U1::T, m0::T, m1::T, λ::T,
        L::_GeometricLagrangian{D, T},
    ) where {D, T}
    h = sqrt(eps(T))
    λp = clamp(λ + h, zero(T), one(T))
    λm = clamp(λ - h, zero(T), one(T))
    return (
        _edge_phi(c_x, c_y0, c_y1, U0, U1, m0, m1, λp, L) -
            _edge_phi(c_x, c_y0, c_y1, U0, U1, m0, m1, λm, L)
    ) / (λp - λm)
end

"""
Edge candidate minimum (§4.4). Evaluates Φ(0), Φ(1), and an ITP root of
dΦ/dλ; the smallest wins. `λ* ∈ {0, 1}` indicates a boundary win.
"""
struct _EdgeDPhi{D, T, LT}
    c_x::SVector{D, T}
    c_y0::SVector{D, T}
    c_y1::SVector{D, T}
    U0::T
    U1::T
    m0::T
    m1::T
    L::LT
end
@inline (g::_EdgeDPhi)(λ) =
    _edge_dphi(g.c_x, g.c_y0, g.c_y1, g.U0, g.U1, g.m0, g.m1, λ, g.L)

@inline function _edge_minimum(
        c_x::SVector{D, T},
        c_y0::SVector{D, T}, c_y1::SVector{D, T},
        U0::T, U1::T, m0::T, m1::T,
        L::_GeometricLagrangian{D, T},
        cutoff::T = T(Inf),
        Φ0::T = _edge_phi(c_x, c_y0, c_y1, U0, U1, m0, m1, zero(T), L),
        Φ1::T = _edge_phi(c_x, c_y0, c_y1, U0, U1, m0, m1, one(T), L),
    ) where {D, T}
    best, λstar = (Φ0 <= Φ1) ? (Φ0, zero(T)) : (Φ1, one(T))
    # Skip ITP if both endpoints already exceed the caller's running best.
    # The interior could still be lower, but on Maupertuis-like problems the
    # interior minimum is bracketed by smaller endpoints in practice.
    (Φ0 > cutoff && Φ1 > cutoff) && return (best, λstar)
    g = _EdgeDPhi(c_x, c_y0, c_y1, U0, U1, m0, m1, L)
    λroot, ok = itp_root(g, zero(T), one(T); atol = T(1.0e-4), maxiter = 20)
    if ok && zero(T) < λroot < one(T)
        Φi = _edge_phi(c_x, c_y0, c_y1, U0, U1, m0, m1, λroot, L)
        if Φi < best
            best, λstar = Φi, λroot
        end
    end
    return (best, λstar)
end

@inline function _edge_slope_estimates(
        state::_OLIMState{D, T},
        y0::CartesianIndex{D}, y1::CartesianIndex{D},
    ) where {D, T}
    dir = y1 - y0
    nbox = state.nbox
    y0m = _shift_in_bounds(y0, -dir, nbox)
    y1p = _shift_in_bounds(y1, dir, nbox)
    m0 = (y0m !== nothing && isfinite(state.U[y0m])) ?
        (state.U[y1] - state.U[y0m]) / 2 : T(NaN)
    m1 = (y1p !== nothing && isfinite(state.U[y1p])) ?
        (state.U[y1p] - state.U[y0]) / 2 : T(NaN)
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
        L::_GeometricLagrangian{3, T}, λ1::T, λ2::T,
    ) where {T}
    y = (one(T) - λ1 - λ2) * c_y0 + λ1 * c_y1 + λ2 * c_y2
    return ((one(T) - λ1 - λ2) * U0 + λ1 * U1 + λ2 * U2) +
        _line_integral(L, y, c_x - y)
end

@inline function _triangle_minimum(
        c_x::SVector{3, T},
        c_y0::SVector{3, T}, c_y1::SVector{3, T}, c_y2::SVector{3, T},
        U0::T, U1::T, U2::T,
        L::_GeometricLagrangian{3, T},
    ) where {T}
    Φ00 = _tri_phi(c_x, c_y0, c_y1, c_y2, U0, U1, U2, L, zero(T), zero(T))
    best, bλ1, bλ2 = Φ00, zero(T), zero(T)
    Φ1 = _tri_phi(c_x, c_y0, c_y1, c_y2, U0, U1, U2, L, one(T), zero(T))
    Φ1 < best && ((best, bλ1, bλ2) = (Φ1, one(T), zero(T)))
    Φ2 = _tri_phi(c_x, c_y0, c_y1, c_y2, U0, U1, U2, L, zero(T), one(T))
    Φ2 < best && ((best, bλ1, bλ2) = (Φ2, zero(T), one(T)))

    Φe01, λs01 = _edge_minimum(c_x, c_y0, c_y1, U0, U1, T(NaN), T(NaN), L)
    if Φe01 < best
        best = Φe01; bλ1 = λs01; bλ2 = zero(T)
    end
    Φe02, λs02 = _edge_minimum(c_x, c_y0, c_y2, U0, U2, T(NaN), T(NaN), L)
    if Φe02 < best
        best = Φe02; bλ1 = zero(T); bλ2 = λs02
    end
    Φe12, λs12 = _edge_minimum(c_x, c_y1, c_y2, U1, U2, T(NaN), T(NaN), L)
    if Φe12 < best
        best = Φe12; bλ1 = one(T) - λs12; bλ2 = λs12
    end

    h = sqrt(eps(T))
    λ1 = one(T) / 3; λ2 = one(T) / 3
    @inline P(a, b) = _tri_phi(c_x, c_y0, c_y1, c_y2, U0, U1, U2, L, a, b)
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
        Φi = _tri_phi(c_x, c_y0, c_y1, c_y2, U0, U1, U2, L, λ1, λ2)
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
        Φ = _vertex_candidate(c_x, c_y, U_y, L)
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
    offsets = _stencil_offsets(Val(K), Val(2))
    cheb = _chebyshev_neighbors(Val(2))
    @inbounds for δ0 in offsets
        y0 = _shift_in_bounds(x, δ0, nbox)
        y0 === nothing && continue
        front[y0] || continue
        U0 = state.U[y0]
        isfinite(U0) || continue
        U0 < best || continue  # vertex bound Φ ≥ U0
        c_y0 = cell_center(grid, y0)
        # Φ0 = vertex candidate from y0 (same as edge Φ(λ=0)). Compute once,
        # reuse for the vertex update and as the Φ(0) endpoint of every edge
        # leaving y0.
        Φ0 = U0 + _line_integral(L, c_y0, c_x - c_y0)
        if Φ0 < best
            best = Φ0
            best_ref = BackRef{2}(y0, y0, NaN32)
        end
        for δch in cheb
            y1 = _shift_in_bounds(y0, δch, nbox)
            y1 === nothing && continue
            y1 > y0 || continue
            front[y1] || continue
            _in_circular_stencil(x, y1, Val(K)) || continue
            U1 = state.U[y1]
            (isfinite(U1) && U1 < best) || continue
            c_y1 = cell_center(grid, y1)
            # Φ(1) is the vertex candidate from y1; cheap, gates ITP.
            Φ1 = U1 + _line_integral(L, c_y1, c_x - c_y1)
            (min(Φ0, Φ1) < best) || continue
            m0, m1 = _edge_slope_estimates(state, y0, y1)
            Φ, λstar = _edge_minimum(
                c_x, c_y0, c_y1, U0, U1, m0, m1, L, best, Φ0, Φ1,
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
    offsets = _stencil_offsets(Val(K), Val(3))
    cheb = _chebyshev_neighbors(Val(3))
    @inbounds for δ0 in offsets
        y0 = _shift_in_bounds(x, δ0, nbox)
        y0 === nothing && continue
        front[y0] || continue
        U0 = state.U[y0]
        isfinite(U0) || continue
        c_y0 = cell_center(grid, y0)
        Φ0 = U0 + _line_integral(L, c_y0, c_x - c_y0)
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
                Φ, λ1, _ = _triangle_minimum(c_x, c_y0, c_y1, c_y2, U0, U1, U2, L)
                if Φ < best
                    best = Φ
                    best_ref = BackRef{3}(y0, y1, Float32(λ1))
                end
            end
        end
    end
    return (best, best_ref)
end
