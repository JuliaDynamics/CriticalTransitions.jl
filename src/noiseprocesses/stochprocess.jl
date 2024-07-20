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
    noise = referrenced_sciml_prob(sys).noise
    covariance = noise.covariance
    if isnothing(covariance)
        return Diagonal(ones(dimension(sys)))
    else
        return covariance
    end
end
