diagonal_cov(l) = Diagonal(ones(l))

"""
$(TYPEDSIGNATURES)

Translates the stochastic process specified in `sys` into the language required by the
`SDEProblem` of `DynamicalSystems.jl`.
"""
function noise_process(sys::CoupledSDEs)
    return sys.integ.W
end

"""
$(TYPEDSIGNATURES)

Gives the covariance matrix specified in `sys`.
"""
function covariance_matrix(sys::CoupledSDEs)
    noise = noise_process(sys)
    covariance = isnothing(noise) ? nothing : noise.covariance
    if isnothing(noise) || isnothing(covariance)
        return diagonal_cov(dimension(sys))
    else
        return covariance
    end
end

"""
$(TYPEDSIGNATURES)

Gives the noise strength specified in `sys`.
"""
function noise_strength(sys::CoupledSDEs)
    sys.noise_strength
end

"""
$(TYPEDSIGNATURES)

Returns the drift field ``b(x)`` of the CoupledSDEs `sys` at the state vector `x`.
"""
function drift(sys::CoupledSDEs{IIP}, x) where {IIP}
    # assumes the drift is time independent
    f = dynamic_rule(sys)
    if IIP
        dx = similar(x)
        f(dx, x, sys.p0, 0)
        return dx
    else
        return f(x, sys.p0, 0)
    end
end
