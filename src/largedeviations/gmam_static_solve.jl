# Preserve static covariance structure in the direct-gMAM metric solve instead of
# materializing a heap Matrix and factorizing it at every interior path point.
function _lambda_theta!(
        θ_out, a_i::StaticArrays.StaticMatrix, b_i, φp, Ainv_b, Ainv_φp, F_cached,
    )
    if F_cached === nothing
        if _isdiag_numerical(a_i)
            @inbounds for k in eachindex(Ainv_b)
                d = a_i[k, k]
                Ainv_b[k] = b_i[k] / d
                Ainv_φp[k] = φp[k] / d
            end
        else
            a_inv = inv(a_i)
            LinearAlgebra.mul!(Ainv_b, a_inv, b_i)
            LinearAlgebra.mul!(Ainv_φp, a_inv, φp)
        end
    else
        copyto!(Ainv_b, b_i)
        LinearAlgebra.ldiv!(F_cached, Ainv_b)
        copyto!(Ainv_φp, φp)
        LinearAlgebra.ldiv!(F_cached, Ainv_φp)
    end

    num = dot(b_i, Ainv_b)
    den = dot(φp, Ainv_φp)
    λ = den > 1.0e-28 ? sqrt(num / den) : zero(eltype(b_i))
    @inbounds for k in eachindex(θ_out)
        θ_out[k] = λ * Ainv_φp[k] - Ainv_b[k]
    end
    return λ
end
