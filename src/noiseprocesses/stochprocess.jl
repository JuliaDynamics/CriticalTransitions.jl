"""
$(TYPEDSIGNATURES)

Translates the stochastic process specified in `sys` into the language required by the
`SDEProblem` of `DynamicalSystems.jl`.
"""
function stochprocess(sys::CoupledSDEs)
    if isnothing(sys.integ.noise)
        prototype = sys.integ.noise_rate_prototype
        if isnothing(prototype) || prototype isa Vector
            return WienerProcess(0.0, zeros(dimension(sys)))
        else
            ArgumentError("Correlated Wiener process not yet implemented.")
            # return CorrelatedWienerProcess(prototype, 0.0, zeros(dimension(sys)))
        end
    else
        return sys.integ.noise
    end
end

"""
$(TYPEDSIGNATURES)

Gives the covariance matrix specified in `sys`.
"""
function covariance_matrix(sys::CoupledSDEs)
    prototype = referrenced_sciml_prob(sys).noise_rate_prototype
    if isnothing(prototype) || prototype isa Vector
        return Diagonal(ones(dimension(sys)))
    else
        ArgumentError("Correlated Wiener process not yet implemented.")
        # return CorrelatedWienerProcess(prototype, 0.0, zeros(dimension(sys)))
    end
end
