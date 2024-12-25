StochasticSystemsBase = Base.get_extension(DynamicalSystemsBase, :StochasticSystemsBase)
covariance_matrix = StochasticSystemsBase.covariance_matrix
diffusion_matrix = StochasticSystemsBase.diffusion_matrix

"""
    StochSystem

An alias to [`CoupledSDEs`](@ref).
This was the name these systems had in CriticalTransitions.jl before v0.3
"""
const StochSystem = CoupledSDEs

"""
$(TYPEDSIGNATURES)

Fetches the stochastic process ``\\mathcal{N}`` specified in the intergrator of `sys`. Returns the type `DiffEqNoiseProcess.NoiseProcess`.
"""
function noise_process(sys::CoupledSDEs)
    return sys.integ.W
end

"""
$(TYPEDSIGNATURES)

Returns the deterministic drift ``f(x)`` of the CoupledSDEs `sys` at state `x`. For time-dependent systems, the time can be specified as a keyword argument `t` (by default `t=0`).
"""
function drift(sys::CoupledSDEs{IIP}, x; t=0) where {IIP}
    f = dynamic_rule(sys)
    if IIP
        dx = similar(x)
        f(dx, x, sys.p0, t)
        return dx
    else
        return f(x, sys.p0, t)
    end
end

"""
$(TYPEDSIGNATURES)

Computes the divergence of the drift field ``f(x)`` at state `x`. For time-
dependent systems, the time can be specified as a keyword argument `t` (by default `t=0`).
"""
function div_drift(sys::CoupledSDEs, x, t=0)
    return tr(jacobian(sys)(x, sys.p0, t))
end;
