function update!(x, xdot, xdotdot, p, pdot, lambda, sys, ϵ; cache = nothing)
    # xdot, p, and lambda are assumed to be pre-computed and consistent with x
    # (via _sgmam_refresh! or central_diff! + update_p!)
    central_diff!(pdot, p)
    Hx = sys.H_x(x, p)
    central_diff!(xdotdot, xdot)
    return update_x!(x, lambda, pdot, xdotdot, Hx, sys, ϵ; cache)
end

function update_x!(x, λ, p′, x′′, Hx, sys::FreidlinWentzellHamiltonian, ϵ; cache = nothing)
    _update_x!(x, λ, p′, x′′, Hx, sys, ϵ, cache)
    return nothing
end

_sgmam_nthreads() =
    hasmethod(Threads.nthreads, Tuple{Symbol}) ?
    (Threads.nthreads(:default) + Threads.nthreads(:interactive)) :
    Threads.nthreads()

# Allocate `nthreads` LinearSolve caches with placeholder Tridiagonal storage of
# the right size. Each `_update_x!` call reuses these by mutating the aliased A
# and b in place and calling `LinearSolve.reinit!`.
function _sgmam_init_tridiag_caches(::Type{T}, L::Int, nthreads::Int) where {T}
    return map(1:nthreads) do _
        d0 = ones(T, L)
        du0 = zeros(T, L - 1)
        dl0 = zeros(T, L - 1)
        T0 = LinearAlgebra.Tridiagonal(dl0, d0, du0)
        rhs0 = zeros(T, L)
        init(
            LinearProblem(T0, rhs0), LUFactorization();
            alias = SciMLBase.LinearAliasSpecifier(; alias_A = true, alias_b = true),
        )
    end
end

# For both Additive and Diagonal noise the implicit step factors per DOF into a
# tridiagonal solve `(I - ϵ λ_i² a_dof,i^{-1} ∂_s²) x_new[dof] = R[dof]`. The only
# difference between the two is whether the coefficient `inv_a(dof, i)` varies in `i`.
# We pass it in as a closure to keep the kernel a single function. `caches` is the
# per-thread LinearSolve cache vector held on the cross-iteration cache struct.
function _update_x_tridiag!(x, λ, p′, x′′, Hx, ϵ, caches, inv_a)
    Nx, Nt = size(x); xa = view(x, :, 1); xb = view(x, :, Nt); idxc = 2:(Nt - 1)
    L = Nt - 2

    Threads.@threads :static for dof in 1:Nx
        cache = caches[Threads.threadid()]
        Tmat = cache.A
        @inbounds for i in 1:L
            ia = inv_a(dof, i + 1)
            Tmat.d[i] = 1 + 2 * ϵ * λ[i + 1]^2 * ia
        end
        @inbounds for i in 1:(L - 1)
            Tmat.du[i] = -ϵ * λ[i + 1]^2 * inv_a(dof, i + 1)
            Tmat.dl[i] = -ϵ * λ[i + 2]^2 * inv_a(dof, i + 2)
        end
        rhs = cache.b
        @inbounds for k in 1:L
            i = k + 1
            ia = inv_a(dof, i)
            rhs[k] = x[dof, i] + ϵ * (λ[i] * p′[dof, i] + Hx[dof, i] - ia * λ[i]^2 * x′′[dof, i])
        end
        rhs[1]   += ϵ * inv_a(dof, 2) * λ[2]^2 * xa[dof]
        rhs[end] += ϵ * inv_a(dof, Nt - 1) * λ[end - 1]^2 * xb[dof]
        LinearSolve.reinit!(cache; A = Tmat, b = rhs)
        solve!(cache)
        x[dof, idxc] .= cache.u
    end
    return nothing
end

"""
Cross-iteration cache for the AdditiveNoise sgMAM `_update_x!` path: holds per-thread
LinearSolve caches and `a_inv` precomputed once (additive `a(x)` is constant).
"""
struct _SgMAMAdditiveCache{T, LC}
    caches::Vector{LC}
    a_inv::Vector{T}
end

"""
Cross-iteration cache for the DiagonalNoise sgMAM `_update_x!` path: holds per-thread
LinearSolve caches and a reusable `a_at` buffer (refilled each call from the current
path snapshot).
"""
struct _SgMAMDiagonalCache{T, LC}
    caches::Vector{LC}
    a_at::Matrix{T}
end

function _update_x!(
        x, λ, p′, x′′, Hx,
        sys::FreidlinWentzellHamiltonian{IIP, D, Hx_t, Hp_t, AF, AdditiveNoise}, ϵ,
        cache::_SgMAMAdditiveCache,
    ) where {IIP, D, Hx_t, Hp_t, AF}
    a_inv = cache.a_inv
    return _update_x_tridiag!(
        x, λ, p′, x′′, Hx, ϵ, cache.caches, (dof, _i) -> a_inv[dof],
    )
end

function _update_x!(
        x, λ, p′, x′′, Hx,
        sys::FreidlinWentzellHamiltonian{IIP, D, Hx_t, Hp_t, AF, DiagonalNoise}, ϵ,
        cache::_SgMAMDiagonalCache,
    ) where {IIP, D, Hx_t, Hp_t, AF}
    Nt = size(x, 2)
    a_at = cache.a_at
    for t in 1:Nt
        a_at[:, t] .= LinearAlgebra.diag(sys.a(view(x, :, t)))
    end
    return _update_x_tridiag!(
        x, λ, p′, x′′, Hx, ϵ, cache.caches, (dof, i) -> 1 / a_at[dof, i],
    )
end

# Fallbacks for direct `update_x!` calls (tests): build a one-shot cache and use it.
function _update_x!(
        x, λ, p′, x′′, Hx,
        sys::FreidlinWentzellHamiltonian{IIP, D, Hx_t, Hp_t, AF, AdditiveNoise}, ϵ,
        ::Nothing,
    ) where {IIP, D, Hx_t, Hp_t, AF}
    Nx, Nt = size(x)
    cache = _build_sgmam_cache(sys, eltype(x), Nx, Nt)
    return _update_x!(x, λ, p′, x′′, Hx, sys, ϵ, cache)
end

function _update_x!(
        x, λ, p′, x′′, Hx,
        sys::FreidlinWentzellHamiltonian{IIP, D, Hx_t, Hp_t, AF, DiagonalNoise}, ϵ,
        ::Nothing,
    ) where {IIP, D, Hx_t, Hp_t, AF}
    Nx, Nt = size(x)
    cache = _build_sgmam_cache(sys, eltype(x), Nx, Nt)
    return _update_x!(x, λ, p′, x′′, Hx, sys, ϵ, cache)
end

"""
Cross-iteration cache for the GeneralNoise sgMAM `_update_x!` path: holds the sparse
block-tridiagonal matrix `M`, an index map from logical (i_in, k1, k2) entries into
`M.nzval`, an RHS buffer, a reusable `LinearSolve.LinearCache`, and an `a_at` buffer
that `_update_p!` fills per iteration and `_update_x!` reads back. The sparsity
pattern is fixed (depends only on `Nt`, `Nx`); only `nzval`, `rhs`, and `a_at` change
per call.
"""
struct _SgMAMGeneralCache{T, LC}
    M::SparseArrays.SparseMatrixCSC{T, Int}
    diag_nzval_idx::Array{Int, 3}    # (Nx, Nx, N_in)
    left_nzval_idx::Matrix{Int}      # (Nx, N_in - 1): for i_in in 2:N_in
    right_nzval_idx::Matrix{Int}     # (Nx, N_in - 1): for i_in in 1:N_in-1
    rhs::Vector{T}
    linear_cache::LC
    a_at::Vector{Matrix{T}}          # length Nt; filled by `_update_p!`
end

"""
Dispatcher that builds a cross-iteration cache for the inner `_update_x!` linear
solve. Each `NoiseShape` gets its own concrete cache so that the `_update_x!`
implementation can avoid rebuilding LinearSolve caches and per-DOF/per-`t` arrays on
every outer iteration.
"""
function _build_sgmam_cache(
        sys::FreidlinWentzellHamiltonian{IIP, D, Hx, Hp, AF, AdditiveNoise},
        ::Type{T}, Nx::Int, Nt::Int,
    ) where {IIP, D, Hx, Hp, AF, T}
    nthreads = _sgmam_nthreads()
    L = Nt - 2
    caches = _sgmam_init_tridiag_caches(T, L, nthreads)
    a_const = sys.a(zeros(T, Nx))
    a_inv = T.(inv.(LinearAlgebra.diag(a_const)))
    return _SgMAMAdditiveCache{T, eltype(caches)}(caches, a_inv)
end

function _build_sgmam_cache(
        ::FreidlinWentzellHamiltonian{IIP, D, Hx, Hp, AF, DiagonalNoise},
        ::Type{T}, Nx::Int, Nt::Int,
    ) where {IIP, D, Hx, Hp, AF, T}
    nthreads = _sgmam_nthreads()
    L = Nt - 2
    caches = _sgmam_init_tridiag_caches(T, L, nthreads)
    a_at = Matrix{T}(undef, Nx, Nt)
    return _SgMAMDiagonalCache{T, eltype(caches)}(caches, a_at)
end

function _build_sgmam_cache(
        ::FreidlinWentzellHamiltonian{IIP, D, Hx, Hp, AF, GeneralNoise},
        ::Type{T}, Nx::Int, Nt::Int,
    ) where {IIP, D, Hx, Hp, AF, T}
    return _build_sgmam_general_cache(T, Nx, Nt)
end

function _find_nzval_idx(M::SparseArrays.SparseMatrixCSC, r::Int, c::Int)
    @inbounds for k in M.colptr[c]:(M.colptr[c + 1] - 1)
        if M.rowval[k] == r
            return k
        end
    end
    error("entry ($r, $c) is not in the sparsity pattern")
end

function _build_sgmam_general_cache(::Type{T}, Nx::Int, Nt::Int) where {T}
    N_in = Nt - 2
    n = N_in * Nx
    nnz_est = N_in * Nx * Nx + 2 * (N_in - 1) * Nx

    Iv = Vector{Int}(undef, nnz_est)
    Jv = Vector{Int}(undef, nnz_est)
    Vv = ones(T, nnz_est)   # placeholder so M is non-singular at build time
    pos = 0
    @inbounds for i_in in 1:N_in
        rb = (i_in - 1) * Nx
        for k1 in 1:Nx, k2 in 1:Nx
            pos += 1; Iv[pos] = rb + k1; Jv[pos] = rb + k2
        end
        if i_in > 1
            base_L = (i_in - 2) * Nx
            for k in 1:Nx
                pos += 1; Iv[pos] = rb + k; Jv[pos] = base_L + k
            end
        end
        if i_in < N_in
            base_R = i_in * Nx
            for k in 1:Nx
                pos += 1; Iv[pos] = rb + k; Jv[pos] = base_R + k
            end
        end
    end
    M = SparseArrays.sparse(Iv, Jv, Vv, n, n)

    diag_idx = Array{Int}(undef, Nx, Nx, N_in)
    left_idx = Matrix{Int}(undef, Nx, N_in - 1)
    right_idx = Matrix{Int}(undef, Nx, N_in - 1)
    @inbounds for i_in in 1:N_in
        rb = (i_in - 1) * Nx
        for k1 in 1:Nx, k2 in 1:Nx
            diag_idx[k1, k2, i_in] = _find_nzval_idx(M, rb + k1, rb + k2)
        end
        if i_in > 1
            base_L = (i_in - 2) * Nx
            for k in 1:Nx
                left_idx[k, i_in - 1] = _find_nzval_idx(M, rb + k, base_L + k)
            end
        end
        if i_in < N_in
            base_R = i_in * Nx
            for k in 1:Nx
                right_idx[k, i_in] = _find_nzval_idx(M, rb + k, base_R + k)
            end
        end
    end

    rhs = zeros(T, n)
    fill!(M.nzval, zero(T))
    lc = init(LinearProblem(M, rhs), LUFactorization())
    a_at = [Matrix{T}(undef, Nx, Nx) for _ in 1:Nt]
    return _SgMAMGeneralCache{T, typeof(lc)}(
        M, diag_idx, left_idx, right_idx, rhs, lc, a_at,
    )
end

# Fallback for direct `update_x!` calls (tests): build a one-shot cache, fill `a_at`,
# then dispatch to the cached variant.
function _update_x!(
        x, λ, p′, x′′, Hx,
        sys::FreidlinWentzellHamiltonian{IIP, D, Hx_t, Hp_t, AF, GeneralNoise}, ϵ,
        ::Nothing = nothing,
    ) where {IIP, D, Hx_t, Hp_t, AF}
    Nx, Nt = size(x)
    cache = _build_sgmam_general_cache(eltype(x), Nx, Nt)
    @inbounds for t in 1:Nt
        a_t = sys.a(view(x, :, t))
        copyto!(cache.a_at[t], a_t isa Matrix ? a_t : Matrix(a_t))
    end
    return _update_x!(x, λ, p′, x′′, Hx, sys, ϵ, cache)
end

function _update_x!(
        x, λ, p′, x′′, Hx,
        sys::FreidlinWentzellHamiltonian{IIP, D, Hx_t, Hp_t, AF, GeneralNoise}, ϵ,
        cache::_SgMAMGeneralCache,
    ) where {IIP, D, Hx_t, Hp_t, AF}
    Nx, Nt = size(x)
    N_in = Nt - 2

    M = cache.M
    rhs = cache.rhs
    xa = view(x, :, 1); xb = view(x, :, Nt)

    @inbounds for i_in in 1:N_in
        i = i_in + 1
        rb = (i_in - 1) * Nx
        # `cache.a_at[i]` was filled by the preceding `_update_p!` call on the same
        # path snapshot, so we skip the closure call here.
        A_i = cache.a_at[i]
        λi2 = λ[i]^2

        for k1 in 1:Nx, k2 in 1:Nx
            v = A_i[k1, k2] + (k1 == k2 ? 2 * ϵ * λi2 : zero(eltype(x)))
            M.nzval[cache.diag_nzval_idx[k1, k2, i_in]] = v
        end
        if i_in > 1
            for k in 1:Nx
                M.nzval[cache.left_nzval_idx[k, i_in - 1]] = -ϵ * λi2
            end
        end
        if i_in < N_in
            for k in 1:Nx
                M.nzval[cache.right_nzval_idx[k, i_in]] = -ϵ * λi2
            end
        end

        for k in 1:Nx
            acc = zero(eltype(x))
            for k2 in 1:Nx
                acc += A_i[k, k2] * (x[k2, i] + ϵ * (λ[i] * p′[k2, i] + Hx[k2, i]))
            end
            v = acc - ϵ * λi2 * x′′[k, i]
            if i_in == 1
                v += ϵ * λi2 * xa[k]
            end
            if i_in == N_in
                v += ϵ * λi2 * xb[k]
            end
            rhs[rb + k] = v
        end
    end

    lc = cache.linear_cache
    LinearSolve.reinit!(lc; A = M, b = rhs)
    solve!(lc)
    sol = lc.u
    @inbounds for i_in in 1:N_in
        rb = (i_in - 1) * Nx
        for k in 1:Nx
            x[k, i_in + 1] = sol[rb + k]
        end
    end
    return nothing
end

function update_p!(p, lambda, x, xdot, sys::FreidlinWentzellHamiltonian; cache = nothing)
    _update_p!(p, lambda, x, xdot, sys, cache)
    return nothing
end

function _update_p!(
        p, lambda, x, xdot,
        sys::FreidlinWentzellHamiltonian{IIP, D, Hx, Hp, AF, AdditiveNoise},
        _cache,
    ) where {IIP, D, Hx, Hp, AF}
    b_ = sys.H_p(x, zero(x))
    a_diag = LinearAlgebra.diag(sys.a(view(x, :, 1)))
    a_inv = inv.(a_diag)
    num = sum(b_ .^ 2 .* a_inv; dims = 1)
    den = sum(xdot .^ 2 .* a_inv; dims = 1)
    lambda .= sqrt.(num ./ den)
    lambda[1] = 0
    lambda[end] = 0
    p .= (lambda .* xdot .- b_) .* a_inv
    return nothing
end

function _update_p!(
        p, lambda, x, xdot,
        sys::FreidlinWentzellHamiltonian{IIP, D, Hx, Hp, AF, DiagonalNoise},
        _cache,
    ) where {IIP, D, Hx, Hp, AF}
    b_ = sys.H_p(x, zero(x))
    Nt = size(x, 2)
    for t in 1:Nt
        a_t = sys.a(view(x, :, t))
        diag_t = LinearAlgebra.diag(a_t)
        inv_t = inv.(diag_t)
        num_t = sum(b_[:, t] .^ 2 .* inv_t)
        den_t = sum(xdot[:, t] .^ 2 .* inv_t)
        λ_t = den_t > 1.0e-28 ? sqrt(num_t / den_t) : zero(eltype(b_))
        lambda[1, t] = isfinite(λ_t) ? λ_t : zero(eltype(b_))
        @inbounds for i in 1:D
            p[i, t] = (lambda[1, t] * xdot[i, t] - b_[i, t]) * inv_t[i]
        end
    end
    lambda[1, 1] = 0
    lambda[1, end] = 0
    return nothing
end

function _update_p!(
        p, lambda, x, xdot,
        sys::FreidlinWentzellHamiltonian{IIP, D, Hx, Hp, AF, GeneralNoise},
        cache,
    ) where {IIP, D, Hx, Hp, AF}
    b_ = sys.H_p(x, zero(x))
    Nt = size(x, 2)
    T = eltype(b_)
    inv_b = Vector{T}(undef, D)
    inv_xdot = Vector{T}(undef, D)
    rhs = Vector{T}(undef, D)
    for t in 1:Nt
        # Evaluate `a_t` once and stash into the GeneralNoise cache so the
        # following `_update_x!` call can skip the closure evaluation.
        a_t = sys.a(view(x, :, t))
        a_t_mat = a_t isa Matrix ? a_t : Matrix(a_t)
        if cache isa _SgMAMGeneralCache
            copyto!(cache.a_at[t], a_t_mat)
        end
        # One LU factorization shared across the three a_t-solves below.
        F = LinearAlgebra.lu(a_t_mat)
        copyto!(inv_b, view(b_, :, t));     LinearAlgebra.ldiv!(F, inv_b)
        copyto!(inv_xdot, view(xdot, :, t)); LinearAlgebra.ldiv!(F, inv_xdot)
        num_t = dot(view(b_, :, t), inv_b)
        den_t = dot(view(xdot, :, t), inv_xdot)
        λ_t = den_t > 1.0e-28 ? sqrt(num_t / den_t) : zero(T)
        lambda[1, t] = λ_t
        @inbounds for i in 1:D
            rhs[i] = λ_t * xdot[i, t] - b_[i, t]
        end
        LinearAlgebra.ldiv!(F, rhs)
        @inbounds for i in 1:D
            p[i, t] = rhs[i]
        end
    end
    lambda[1, 1] = 0
    lambda[1, end] = 0
    return nothing
end

function central_diff!(xdot, x)
    # ̇xₙ = 0.5(xₙ₊₁ - xₙ₋₁) central finite difference
    @views xdot[:, 2:(end - 1)] .= 0.5 .* (x[:, 3:end] .- x[:, 1:(end - 2)])
    return nothing
end

function _sgmam_refresh!(xdot, p, lambda, x, sys; cache = nothing)
    central_diff!(xdot, x)
    update_p!(p, lambda, x, xdot, sys; cache)
    return nothing
end

# Freidlin-Wentzell action of an instanton path discretized over a uniform arclength
# parameter `s ∈ [0, 1]`. With `xdot` produced by `central_diff!` (which does not
# divide by Δs), `xdot[:, i] ≈ Δs · dϕ/ds`, so `dot(xdot, p) ≈ ∫ p · dϕ`. On the
# zero-energy shell `H = p·f + |p|²/2 = 0`, the identity `p · dϕ = (|p|²/2) dt`
# holds along the instanton, so `∫ p · dϕ = ½ ∫ |p|² dt = S_FW`. No additional `/2`.
FW_action(xdot, p) = dot(xdot, p)
