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

# add covariance to noise if it is nothing (resolves ambiguity)
function resolve_covariance!(prob, noise_type, param, D, IIP)
    noise = prob.noise
    noise_prototype = prob.noise_rate_prototype
    noise_size = isnothing(noise_prototype) ? nothing : size(noise_prototype)
    encoded_cov = noise_type[:additive] && noise_type[:autonomous] && noise_size == (D, D)
    if encoded_cov && !isnothing(noise) && isnothing(noise.covariance)
        if IIP
            du = deepcopy(noise_prototype)
            prob.g(du, zeros(D), param, 0.0)
            noise.covariance = du
        else
            noise.covariance = prob.g(zeros(D), param, 0.0)
        end # NoiseProcess is a mutable struct
    end
end

"""
$(TYPEDSIGNATURES)

Gives the noise strength specified in `sys`.
"""
function noise_strength(sys::CoupledSDEs)
    return sys.noise_strength
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

"""
$(TYPEDSIGNATURES)

Computes the divergence of the drift field `sys.f` at the given point `x`.
"""
function div_drift(sys::CoupledSDEs, x)
    b(x) = drift(sys, x)
    return tr(ForwardDiff.jacobian(b, x))
end;
