struct SgMAMDecoupledCache{T, LC}
    a_inv::Matrix{T}                # diag(a⁻¹) per (k, t)
    Ainv_b::Matrix{T}               # a⁻¹ · b stored elementwise
    Ainv_xd::Matrix{T}              # a⁻¹ · ẋ stored elementwise
    p_zero::Matrix{T}               # zero buffer reused as 2nd arg to sys.H_p
    Hp_buf::Matrix{T}               # buffer for H_p output (Nx × Nt)
    Hx_buf::Matrix{T}               # buffer for H_x output (Nx × Nt)
    linear_cache::LC
end

struct SgMAMCoupledCache{T, LC, LU}
    M::SparseArrays.SparseMatrixCSC{T, Int}
    diag_idx::Array{Int, 3}
    off_idx::Matrix{Int}
    linear_cache::LC
    a_at::Vector{Matrix{T}}         # a(x_t) stored per t (only mutated for state-dep a)
    a_const_lu::LU                  # LU(a) for constant a; `nothing` for state-dep
    Hp_buf::Matrix{T}
    Hx_buf::Matrix{T}
    rhs::Vector{T}
    Ainv_b::Vector{T}               # a(x_t)⁻¹ · b_t  (scratch, per t)
    Ainv_xd::Vector{T}              # a(x_t)⁻¹ · ẋ_t  (scratch, per t)
    rhs_p::Vector{T}                # scratch for p update
    p_zero::Matrix{T}               # zero buffer reused as 2nd arg to sys.H_p
end

function _init_tridiag_cache(::Type{T}, L::Int) where {T}
    dl = zeros(T, L - 1); d = ones(T, L); du = zeros(T, L - 1)
    Tmat = LinearAlgebra.Tridiagonal(dl, d, du)
    rhs = zeros(T, L)
    return init(
        LinearProblem(Tmat, rhs), LUFactorization();
        alias = SciMLBase.LinearAliasSpecifier(; alias_A = true, alias_b = true),
    )
end

function _fill_constant_a_inv!(a_inv::Matrix{T}, a, Nx, Nt) where {T}
    diag_inv = T.(inv.(LinearAlgebra.diag(a(zeros(T, Nx)))))
    @inbounds for t in 1:Nt, k in 1:Nx
        a_inv[k, t] = diag_inv[k]
    end
    return nothing
end

function _refill_state_dep_a_inv!(a_inv::Matrix{T}, a, x) where {T}
    @inbounds for t in axes(x, 2)
        d = LinearAlgebra.diag(a(view(x, :, t)))
        for k in axes(a_inv, 1)
            a_inv[k, t] = inv(d[k])
        end
    end
    return nothing
end

function _build_decoupled_cache(sys, ::Type{T}, Nx::Int, Nt::Int) where {T}
    a_inv = Matrix{T}(undef, Nx, Nt)
    Ainv_b = Matrix{T}(undef, Nx, Nt)
    Ainv_xd = Matrix{T}(undef, Nx, Nt)
    p_zero = zeros(T, Nx, Nt)
    is_constant(sys.a) === Val(true) && _fill_constant_a_inv!(a_inv, sys.a, Nx, Nt)
    linear_cache = _init_tridiag_cache(T, Nt - 2)
    Hp_buf = Matrix{T}(undef, Nx, Nt)
    Hx_buf = Matrix{T}(undef, Nx, Nt)
    return SgMAMDecoupledCache{T, typeof(linear_cache)}(
        a_inv, Ainv_b, Ainv_xd, p_zero, Hp_buf, Hx_buf, linear_cache,
    )
end

function _find_nzval_idx(M::SparseArrays.SparseMatrixCSC, r::Int, c::Int)
    @inbounds for k in M.colptr[c]:(M.colptr[c + 1] - 1)
        if M.rowval[k] == r
            return k
        end
    end
    return error("entry ($r, $c) is not in the sparsity pattern")
end

function _build_coupled_cache(sys, ::Type{T}, Nx::Int, Nt::Int) where {T}
    N_in = Nt - 2
    n = N_in * Nx
    nnz_est = N_in * Nx^2 + 2 * (N_in - 1) * Nx
    Iv = Vector{Int}(undef, nnz_est)
    Jv = Vector{Int}(undef, nnz_est)
    Vv = ones(T, nnz_est)
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

    diag_idx = Array{Int, 3}(undef, Nx, Nx, N_in)
    off_idx = Matrix{Int}(undef, Nx, 2 * (N_in - 1))
    @inbounds for i_in in 1:N_in
        rb = (i_in - 1) * Nx
        for k1 in 1:Nx, k2 in 1:Nx
            diag_idx[k1, k2, i_in] = _find_nzval_idx(M, rb + k1, rb + k2)
        end
        if i_in > 1
            base_L = (i_in - 2) * Nx
            for k in 1:Nx
                off_idx[k, 2 * (i_in - 1) - 1] = _find_nzval_idx(M, rb + k, base_L + k)
            end
        end
        if i_in < N_in
            base_R = i_in * Nx
            for k in 1:Nx
                off_idx[k, 2 * (i_in - 1) + 2] = _find_nzval_idx(M, rb + k, base_R + k)
            end
        end
    end

    rhs = zeros(T, n)
    fill!(M.nzval, zero(T))
    lc = init(
        LinearProblem(M, rhs), LUFactorization();
        alias = SciMLBase.LinearAliasSpecifier(; alias_A = true, alias_b = true),
    )
    a_at = [Matrix{T}(undef, Nx, Nx) for _ in 1:Nt]
    p_zero = zeros(T, Nx, Nt)
    Hp_buf = Matrix{T}(undef, Nx, Nt)
    Hx_buf = Matrix{T}(undef, Nx, Nt)
    # Precompute LU(a) once for constant a (same for all t).
    a_const_lu = if is_constant(sys.a) === Val(true)
        a0 = sys.a(zeros(T, Nx))
        a0M = a0 isa Matrix ? Matrix{T}(a0) : Matrix{T}(a0)
        for t in 1:Nt
            copyto!(a_at[t], a0M)
        end
        LinearAlgebra.lu(a0M)
    else
        nothing
    end
    return SgMAMCoupledCache{T, typeof(lc), typeof(a_const_lu)}(
        M, diag_idx, off_idx, lc, a_at, a_const_lu, Hp_buf, Hx_buf, rhs,
        Vector{T}(undef, Nx), Vector{T}(undef, Nx), Vector{T}(undef, Nx),
        p_zero,
    )
end

function build_sgmam_cache(
        sys::FreidlinWentzellHamiltonian, x_initial::AbstractMatrix{T}, Nt::Int,
    ) where {T}
    Nx = size(x_initial, 1)
    x_ref = collect(view(x_initial, :, 1))
    is_diagonal = _validate_and_classify_a(sys.a, x_ref)
    return _build_sgmam_cache(Val(is_diagonal), sys, T, Nx, Nt)
end

_build_sgmam_cache(::Val{true}, sys, T, Nx, Nt) = _build_decoupled_cache(sys, T, Nx, Nt)
_build_sgmam_cache(::Val{false}, sys, T, Nx, Nt) = _build_coupled_cache(sys, T, Nx, Nt)

function update_p!(p, λ, x, xdot, sys, cache::SgMAMDecoupledCache)
    if is_constant(sys.a) === Val(false)
        _refill_state_dep_a_inv!(cache.a_inv, sys.a, x)
    end
    b_ = _eval_Hp!(cache.Hp_buf, sys, x, cache.p_zero)
    @. cache.Ainv_b = b_ * cache.a_inv
    @. cache.Ainv_xd = xdot * cache.a_inv
    num = sum(b_ .* cache.Ainv_b; dims = 1)
    den = sum(xdot .* cache.Ainv_xd; dims = 1)
    @. λ = ifelse(den > 1.0e-28, sqrt(num / den), zero(eltype(num)))
    λ[1, 1] = λ[1, end] = 0
    @. p = (λ * xdot - b_) * cache.a_inv
    return nothing
end

function update_p!(p, λ, x, xdot, sys, cache::SgMAMCoupledCache)
    b_ = _eval_Hp!(cache.Hp_buf, sys, x, cache.p_zero)
    @inbounds for t in axes(x, 2)
        F = _coupled_lu_at_t!(cache, sys, view(x, :, t), t)
        bt = view(b_, :, t); xt = view(xdot, :, t)
        copyto!(cache.Ainv_b, bt); LinearAlgebra.ldiv!(F, cache.Ainv_b)
        copyto!(cache.Ainv_xd, xt); LinearAlgebra.ldiv!(F, cache.Ainv_xd)
        num = dot(bt, cache.Ainv_b)
        den = dot(xt, cache.Ainv_xd)
        λt = den > 1.0e-28 ? sqrt(num / den) : zero(eltype(b_))
        λ[1, t] = λt
        for k in axes(p, 1)
            cache.rhs_p[k] = λt * xt[k] - bt[k]
        end
        LinearAlgebra.ldiv!(F, cache.rhs_p)
        for k in axes(p, 1)
            p[k, t] = cache.rhs_p[k]
        end
    end
    λ[1, 1] = λ[1, end] = 0
    return nothing
end

# Return LU(a(x_t)). For constant `a` reuse the prebuilt factorization stored in
# `cache.a_const_lu`; otherwise refresh `cache.a_at[t]` from `sys.a(x_t)` and factor it.
@inline _coupled_lu_at_t!(cache, sys, _xt, _t, F::LinearAlgebra.LU) = F
@inline function _coupled_lu_at_t!(cache, sys, xt, t, ::Nothing)
    a_t = sys.a(xt)
    copyto!(cache.a_at[t], a_t isa Matrix ? a_t : Matrix(a_t))
    return LinearAlgebra.lu(cache.a_at[t])
end
@inline _coupled_lu_at_t!(cache, sys, xt, t) =
    _coupled_lu_at_t!(cache, sys, xt, t, cache.a_const_lu)

# update_x! relies on update_p! (via _sgmam_refresh!) having populated cache.a_inv /
# cache.a_at for the current x. No re-refill here.
function update_x!(
        x, λ, p′, x′′, Hx, sys::FreidlinWentzellHamiltonian, ϵ,
        cache::SgMAMDecoupledCache,
    )
    a_inv = cache.a_inv
    lc = cache.linear_cache
    Nx, Nt = size(x); L = Nt - 2
    xa = view(x, :, 1); xb = view(x, :, Nt)
    Tmat = lc.A
    rhs = lc.b
    for dof in 1:Nx
        @inbounds for i in 1:L
            Tmat.d[i] = 1 + 2 * ϵ * λ[i + 1]^2 * a_inv[dof, i + 1]
        end
        @inbounds for i in 1:(L - 1)
            Tmat.du[i] = -ϵ * λ[i + 1]^2 * a_inv[dof, i + 1]
            Tmat.dl[i] = -ϵ * λ[i + 2]^2 * a_inv[dof, i + 2]
        end
        @inbounds for k in 1:L
            i = k + 1
            ia = a_inv[dof, i]
            rhs[k] = x[dof, i] + ϵ * (
                λ[i] * p′[dof, i] + Hx[dof, i] - ia * λ[i]^2 * x′′[dof, i]
            )
        end
        rhs[1] += ϵ * a_inv[dof, 2] * λ[2]^2 * xa[dof]
        rhs[end] += ϵ * a_inv[dof, Nt - 1] * λ[end - 1]^2 * xb[dof]
        LinearSolve.reinit!(lc; A = Tmat, b = rhs)
        solve!(lc)
        @inbounds for k in 1:L
            x[dof, k + 1] = lc.u[k]
        end
    end
    return nothing
end

function update_x!(
        x, λ, p′, x′′, Hx, sys::FreidlinWentzellHamiltonian, ϵ,
        cache::SgMAMCoupledCache,
    )
    Nx, Nt = size(x); N_in = Nt - 2
    M = cache.M
    rhs = cache.rhs
    xa = view(x, :, 1); xb = view(x, :, Nt)
    @inbounds for i_in in 1:N_in
        i = i_in + 1
        rb = (i_in - 1) * Nx
        A_i = cache.a_at[i]
        λi2 = λ[i]^2
        for k1 in 1:Nx, k2 in 1:Nx
            v = A_i[k1, k2] + (k1 == k2 ? 2 * ϵ * λi2 : zero(eltype(x)))
            M.nzval[cache.diag_idx[k1, k2, i_in]] = v
        end
        if i_in > 1
            for k in 1:Nx
                M.nzval[cache.off_idx[k, 2 * (i_in - 1) - 1]] = -ϵ * λi2
            end
        end
        if i_in < N_in
            for k in 1:Nx
                M.nzval[cache.off_idx[k, 2 * (i_in - 1) + 2]] = -ϵ * λi2
            end
        end
        b_vec = cache.Ainv_b
        for k2 in 1:Nx
            b_vec[k2] = x[k2, i] + ϵ * (λ[i] * p′[k2, i] + Hx[k2, i])
        end
        rhs_view = view(rhs, (rb + 1):(rb + Nx))
        LinearAlgebra.mul!(rhs_view, A_i, b_vec)
        for k in 1:Nx
            rhs_view[k] -= ϵ * λi2 * x′′[k, i]
            i_in == 1    && (rhs_view[k] += ϵ * λi2 * xa[k])
            i_in == N_in && (rhs_view[k] += ϵ * λi2 * xb[k])
        end
    end
    LinearSolve.reinit!(cache.linear_cache; A = M, b = rhs)
    solve!(cache.linear_cache)
    sol = cache.linear_cache.u
    @inbounds for i_in in 1:N_in
        rb = (i_in - 1) * Nx
        for k in 1:Nx
            x[k, i_in + 1] = sol[rb + k]
        end
    end
    return nothing
end

function central_diff!(xdot, x)
    @views xdot[:, 2:(end - 1)] .= 0.5 .* (x[:, 3:end] .- x[:, 1:(end - 2)])
    return nothing
end

function _sgmam_refresh!(xdot, p, lambda, x, sys, cache)
    central_diff!(xdot, x)
    update_p!(p, lambda, x, xdot, sys, cache)
    return nothing
end

function update!(x, xdot, xdotdot, p, pdot, lambda, sys, ϵ, cache)
    central_diff!(pdot, p)
    Hx = _eval_Hx!(cache.Hx_buf, sys, x, p)
    central_diff!(xdotdot, xdot)
    return update_x!(x, lambda, pdot, xdotdot, Hx, sys, ϵ, cache)
end

FW_action(xdot, p) = dot(xdot, p)
