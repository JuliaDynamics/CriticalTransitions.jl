StochasticSystemsBase = Base.get_extension(DynamicalSystemsBase, :StochasticSystemsBase)
covariance_matrix = StochasticSystemsBase.covariance_matrix
diffusion_matrix = StochasticSystemsBase.diffusion_matrix

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
function drift(sys::Union{CoupledSDEs{IIP}, CoupledODEs{IIP}}, x; t = 0) where {IIP}
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
function div_drift(sys::ContinuousTimeDynamicalSystem, x, t = 0)
    return tr(jacobian(sys)(x, sys.p0, t))
end;

"""
$(TYPEDSIGNATURES)

Returns the SDE solver specified in the `diffeq` settings of the `CoupledSDEs`.
"""
solver(ds::ContinuousTimeDynamicalSystem) = ds.integ.alg

"""
$(TYPEDSIGNATURES)

Returns the effective noise strength ``\\sigma_{\\mathrm{eff}}`` of the `CoupledSDEs`
`sys`, defined as ``\\sigma_{\\mathrm{eff}} = \\sqrt{\\mathrm{tr}(\\Sigma)/D}`` where
``\\Sigma = `` `covariance_matrix(sys)` is the SDE's noise rate matrix and `D` its
state-space dimension.

This is the σ produced by the trace-normalization convention used in the large-deviation
action functionals (see [`fw_action`](@ref)): the SDE
``\\mathrm{d}x = b\\,\\mathrm{d}t + \\sigma\\,\\Sigma_0\\,\\mathrm{d}W`` is canonicalized as
``\\sigma_{\\mathrm{eff}}^2 = \\mathrm{tr}(\\sigma^2\\,Q)/D`` with
``Q = \\Sigma_0\\Sigma_0^\\top``. For an SDE constructed with isotropic `covariance = I`,
`noise_strength(sys)` recovers the `noise_strength` keyword passed at construction time.
More generally, the returned value is the per-direction average of the noise variance.

For an SDE built with a custom `g` that is not of the form `σ * sqrt(Q)`, the returned
value is still the trace-based effective σ of the resulting covariance matrix.
"""
function noise_strength(sys::CoupledSDEs)
    Σ = covariance_matrix(sys)
    return sqrt(tr(Σ) / size(Σ, 1))
end
