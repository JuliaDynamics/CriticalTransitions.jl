# Preserve static covariance structure in the coupled sgMAM momentum solve instead of
# materializing and factorizing a heap Matrix at every path point.
struct _StaticCoupledCovarianceFactor{A}
    a_inv::A
end

@inline function _mul_static_inverse!(
        b, a_inv::StaticArrays.StaticMatrix{N, N},
    ) where {N}
    rhs = ntuple(i -> @inbounds(b[i]), Val(N))
    @inbounds for i in 1:N
        acc = zero(eltype(b))
        for j in 1:N
            acc += a_inv[i, j] * rhs[j]
        end
        b[i] = acc
    end
    return b
end

@inline function LinearAlgebra.ldiv!(F::_StaticCoupledCovarianceFactor, b::AbstractVector)
    return _mul_static_inverse!(b, F.a_inv)
end

@inline _state_dependent_coupled_factor(a::StaticArrays.StaticMatrix, _a_buf) =
    _StaticCoupledCovarianceFactor(inv(a))

@inline _state_dependent_coupled_factor(_a, a_buf) = LinearAlgebra.lu(a_buf)

# More-specific state-dependent dispatch for coupled sgMAM. `a_at[t]` is still filled
# exactly as before because the sparse x-update consumes the dense covariance buffers.
@inline function _coupled_lu_at_t!(
        cache::SgMAMCoupledCache, sys::FreidlinWentzellHamiltonian, xt, t, ::Nothing,
    )
    a_t = sys.a(xt)
    copyto!(cache.a_at[t], a_t)
    return _state_dependent_coupled_factor(a_t, cache.a_at[t])
end
