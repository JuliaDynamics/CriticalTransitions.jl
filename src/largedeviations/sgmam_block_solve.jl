# Exploit the block-tridiagonal structure of the coupled sgMAM x-update. The existing
# sparse solver remains the fallback for non-symmetric user metrics or a failed block factorization.
struct SgMAMBlockCoupledCache{C, S, F, R, I, V}
    base::C
    schur::S
    factors::F
    rhs_blocks::R
    inv_prev::I
    tmp::V
end

function _build_coupled_cache(
        sys::FreidlinWentzellHamiltonian, ::Type{T}, Nx::Int, Nt::Int,
    ) where {T}
    base = invoke(
        _build_coupled_cache, Tuple{Any, Type{T}, Int, Int}, sys, T, Nx, Nt,
    )
    L = Nt - 2
    schur = [Matrix{T}(LinearAlgebra.I, Nx, Nx) for _ in 1:L]
    factors = [
        LinearAlgebra.cholesky!(LinearAlgebra.Hermitian(S, :L); check = false) for
            S in schur
    ]
    rhs_blocks = Matrix{T}(undef, Nx, L)
    inv_prev = Matrix{T}(undef, Nx, Nx)
    tmp = Vector{T}(undef, Nx)
    return SgMAMBlockCoupledCache(
        base, schur, factors, rhs_blocks, inv_prev, tmp,
    )
end

function update_p!(p, lambda, x, xdot, sys, cache::SgMAMBlockCoupledCache)
    return update_p!(p, lambda, x, xdot, sys, cache.base)
end

function _sgmam_block_applicable(cache::SgMAMBlockCoupledCache)
    @inbounds for i_in in eachindex(cache.schur)
        LinearAlgebra.issymmetric(cache.base.a_at[i_in + 1]) || return false
    end
    return true
end

function _block_thomas_solve!(
        schur, factors, rhs_blocks, inv_prev, tmp, eps_step, lambda,
    )
    Nx, L = size(rhs_blocks)
    eps2 = eps_step * eps_step
    @inbounds for i in 1:L
        S = schur[i]
        if i > 1
            q_i = lambda[i + 1]^2
            q_prev = lambda[i]^2
            Cprev = factors[i - 1]
            fill!(inv_prev, zero(eltype(inv_prev)))
            for k in 1:Nx
                inv_prev[k, k] = one(eltype(inv_prev))
            end
            LinearAlgebra.ldiv!(Cprev, inv_prev)
            coupling2 = eps2 * q_i * q_prev
            for k2 in 1:Nx, k1 in 1:Nx
                S[k1, k2] -= coupling2 * inv_prev[k1, k2]
            end
            for k in 1:Nx
                tmp[k] = rhs_blocks[k, i - 1]
            end
            LinearAlgebra.ldiv!(Cprev, tmp)
            coupling = eps_step * q_i
            for k in 1:Nx
                rhs_blocks[k, i] += coupling * tmp[k]
            end
        end
        C = LinearAlgebra.cholesky!(LinearAlgebra.Hermitian(S, :L); check = false)
        LinearAlgebra.issuccess(C) || return false
    end

    @inbounds for i in L:-1:1
        coupling = eps_step * lambda[i + 1]^2
        for k in 1:Nx
            tmp[k] = rhs_blocks[k, i]
            i < L && (tmp[k] += coupling * rhs_blocks[k, i + 1])
        end
        LinearAlgebra.ldiv!(factors[i], tmp)
        for k in 1:Nx
            rhs_blocks[k, i] = tmp[k]
        end
    end
    return true
end

function update_x!(
        x, lambda, pdot, xdotdot, Hx, sys::FreidlinWentzellHamiltonian, eps_step,
        cache::SgMAMBlockCoupledCache,
    )
    _sgmam_block_applicable(cache) ||
        return update_x!(x, lambda, pdot, xdotdot, Hx, sys, eps_step, cache.base)

    base = cache.base
    Nx, Nt = size(x)
    L = Nt - 2
    xa = view(x, :, 1)
    xb = view(x, :, Nt)
    @inbounds for i_in in 1:L
        i = i_in + 1
        q = lambda[i]^2
        A_i = base.a_at[i]
        S = cache.schur[i_in]
        for k2 in 1:Nx, k1 in 1:Nx
            S[k1, k2] = A_i[k1, k2] +
                (k1 == k2 ? 2 * eps_step * q : zero(eltype(x)))
        end

        b_vec = base.Ainv_b
        for k in 1:Nx
            b_vec[k] = x[k, i] + eps_step * (lambda[i] * pdot[k, i] + Hx[k, i])
        end
        rhs_i = view(cache.rhs_blocks, :, i_in)
        LinearAlgebra.mul!(rhs_i, A_i, b_vec)
        coupling = eps_step * q
        for k in 1:Nx
            rhs_i[k] -= coupling * xdotdot[k, i]
            i_in == 1 && (rhs_i[k] += coupling * xa[k])
            i_in == L && (rhs_i[k] += coupling * xb[k])
        end
    end

    _block_thomas_solve!(
        cache.schur, cache.factors, cache.rhs_blocks, cache.inv_prev, cache.tmp,
        eps_step, lambda,
    ) || return update_x!(x, lambda, pdot, xdotdot, Hx, sys, eps_step, base)

    @inbounds for i_in in 1:L, k in 1:Nx
        x[k, i_in + 1] = cache.rhs_blocks[k, i_in]
    end
    return nothing
end

function update!(
        x, xdot, xdotdot, p, pdot, lambda, sys, eps_step,
        cache::SgMAMBlockCoupledCache,
    )
    central_diff!(pdot, p)
    Hx = _eval_Hx!(cache.base.Hx_buf, sys, x, p)
    central_diff!(xdotdot, xdot)
    return update_x!(x, lambda, pdot, xdotdot, Hx, sys, eps_step, cache)
end
