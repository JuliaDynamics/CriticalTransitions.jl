"""
Solve `A'P + PA + 2 P Q P = 0` for symmetric positive definite `P` via the
symplectic Hamiltonian `H = [-A -2Q; 0 A']` and extract the stable
invariant subspace.
"""
function _care(A::AbstractMatrix{T}, Q::AbstractMatrix{T}) where {T}
    D = LinearAlgebra.checksquare(A)
    H = [
        Matrix(-A)        Matrix(-2 * Q);
        zeros(T, D, D)     Matrix(A')
    ]
    F = eigen(H)
    stable = findall(λ -> real(λ) < -sqrt(eps(real(T))), F.values)
    length(stable) == D || throw(
        ArgumentError(
            "CARE: drift Jacobian has no $D-dim stable invariant subspace",
        ),
    )
    V = F.vectors[:, stable]
    X1 = V[1:D, :]; X2 = V[(D + 1):(2D), :]
    P = real(X2 / X1)
    return (P + P') / 2
end

function _seed_near_source!(
        state::_OLIMState{D, T},
        grid::CartesianGrid{D, T},
        source::CartesianIndex{D},
        sys::CoupledSDEs,
        L::_GeometricLagrangian{D, T},
        ::Val{K_seed};
        verbose::Bool = false,
    ) where {D, T, K_seed}
    state.U[source] = zero(T)
    state.status[source] = _ACCEPTED
    K_seed == 0 && return state

    x_A = cell_center(grid, source)
    A_jac = ForwardDiff.jacobian(L.b, x_A)
    λ = eigen(A_jac).values
    if any(λi -> real(λi) >= -sqrt(eps(T)), λ)
        verbose &&
            @info "OLIM: attractor cell is not strictly stable; analytic seed skipped"
        return state
    end
    Q_mat = _diffusion_at(L, x_A)
    P = try
        _care(Matrix(A_jac), Matrix(Q_mat))
    catch err
        verbose && @info "OLIM: CARE solve failed ($(err)); skipping analytic seed"
        return state
    end

    nbox = state.nbox
    axes_seed = ntuple(_ -> -K_seed:K_seed, Val(D))
    @inbounds for δ in CartesianIndices(axes_seed)
        iszero(δ) && continue
        I = source + δ
        in_bounds = true
        for d in 1:D
            (1 <= I[d] <= nbox[d]) || (in_bounds = false; break)
        end
        in_bounds || continue
        x = cell_center(grid, I)
        dx = x - x_A
        state.U[I] = dot(dx, P * dx)
        state.back_pointer[I] = BackRef{D}(source, source, NaN32)
        state.status[I] = _ACCEPTED
    end
    return state
end

"""
Dijkstra-style one-pass sweep filling `state.U`. After analytic seeding,
flips seed cells that touch UNKNOWN to FRONT and seeds the heap with their
first-tier CONSIDERED neighbours, then iterates until the heap is empty.
"""
function _sweep!(
        state::_OLIMState{D, T},
        grid::CartesianGrid{D, T},
        source::CartesianIndex{D},
        sys::CoupledSDEs,
        L::_GeometricLagrangian{D, T},
        ::Val{K}, ::Val{K_seed};
        verbose::Bool, show_progress::Bool,
    ) where {D, T, K, K_seed}
    _seed_near_source!(state, grid, source, sys, L, Val(K_seed); verbose = verbose)

    nbox = state.nbox
    N = prod(nbox)
    cheb = _chebyshev_neighbors(Val(D))

    @inbounds for x in CartesianIndices(state.U)
        state.status[x] == _ACCEPTED || continue
        touches_unknown = false
        for δ in cheb
            n = _shift_in_bounds(x, δ, nbox); n === nothing && continue
            state.status[n] == _UNKNOWN && (touches_unknown = true; break)
        end
        if touches_unknown
            state.status[x] = _FRONT
            state.front[x] = true
        end
    end

    @inbounds for x in CartesianIndices(state.U)
        state.status[x] == _FRONT || continue
        for δ in _stencil_offsets(Val(K), Val(D))
            n = _shift_in_bounds(x, δ, nbox); n === nothing && continue
            state.status[n] == _UNKNOWN || continue
            Φ, br = _local_update(Val(D), n, state, grid, L, Val(K))
            isfinite(Φ) || continue
            if Φ < state.U[n]
                state.U[n] = Φ
                state.back_pointer[n] = br
            end
            lin = _linear(n, nbox)
            if state.status[n] == _UNKNOWN
                state.status[n] = _CONSIDERED
                state.handles[lin] = push!(state.heap, (state.U[n], lin))
            end
        end
    end

    if show_progress
        return _sweep_loop!(
            state, grid, L, Val(K), Val(D),
            Progress(N; desc = "OLIM sweep")
        )
    else
        return _sweep_loop!(state, grid, L, Val(K), Val(D), nothing)
    end
end

@inline _maybe_tick(::Nothing) = nothing
@inline _maybe_tick(p::Progress) = next!(p)

function _sweep_loop!(
        state::_OLIMState{D, T}, grid::CartesianGrid{D, T},
        L::_GeometricLagrangian{D, T},
        ::Val{K}, ::Val{D}, prog,
    ) where {D, T, K}
    nbox = state.nbox
    cheb = _chebyshev_neighbors(Val(D))
    while !isempty(state.heap)
        (_, lin_c) = pop!(state.heap)
        c = CartesianIndices(nbox)[lin_c]
        state.handles[lin_c] = 0
        state.status[c] = _FRONT
        state.front[c] = true
        _maybe_tick(prog)

        for δ in _stencil_offsets(Val(K), Val(D))
            xnew = _shift_in_bounds(c, δ, nbox); xnew === nothing && continue
            state.status[xnew] == _ACCEPTED && continue
            state.status[xnew] == _FRONT    && continue
            Φ, br = _local_update(Val(D), xnew, state, grid, L, Val(K))
            (isfinite(Φ) && Φ < state.U[xnew]) || continue
            state.U[xnew] = Φ
            state.back_pointer[xnew] = br
            lin_x = _linear(xnew, nbox)
            if state.status[xnew] == _UNKNOWN
                state.status[xnew] = _CONSIDERED
                state.handles[lin_x] = push!(state.heap, (Φ, lin_x))
            else
                DataStructures.update!(state.heap, state.handles[lin_x], (Φ, lin_x))
            end
        end

        prune = true
        for δ in cheb
            n = _shift_in_bounds(c, δ, nbox); n === nothing && continue
            s = state.status[n]
            if s == _UNKNOWN || s == _CONSIDERED
                prune = false; break
            end
        end
        if prune
            state.status[c] = _ACCEPTED
            state.front[c] = false
        end
    end
    return state
end
