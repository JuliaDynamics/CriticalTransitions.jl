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

Effective noise strength ``\\sigma_{\\mathrm{eff}} = \\sqrt{\\mathrm{tr}(\\mathbf{\\Sigma})/D}``
of `sys`, where ``\\mathbf{\\Sigma} = `` `covariance_matrix(sys)` and ``D = `` `dimension(sys)`.
Together with [`normalize_covariance!`](@ref) it satisfies
``\\sigma_{\\mathrm{eff}}^2 \\cdot \\mathbf{Q}_{\\mathrm{can}} = \\mathbf{\\Sigma}``.

For an SDE built as ``\\mathrm{d}\\mathbf{x} = \\mathbf{b}\\,\\mathrm{d}t + \\sigma\\sqrt{\\mathbf{Q}}\\,\\mathrm{d}\\mathbf{W}``:

  - **Isotropic `covariance = c·I`**: ``\\sigma_{\\mathrm{eff}} = \\sigma\\sqrt{c}``. The
    default `covariance = I` (``c=1``) recovers the construction-time `noise_strength`
    keyword exactly.
  - **Anisotropic `Q`**: ``\\sigma_{\\mathrm{eff}} = \\sigma\\sqrt{\\mathrm{tr}(\\mathbf{Q})/D}``,
    the per-direction average noise amplitude. Equals ``\\sigma`` whenever
    ``\\mathrm{tr}(\\mathbf{Q}) = D``.
  - **Custom `g`** (not of the form ``\\sigma\\sqrt{\\mathbf{Q}}``): still well-defined as the
    trace-based effective σ of the resulting ``\\mathbf{\\Sigma}``.

See [Large deviation theory](@ref) for the role of ``\\sigma_{\\mathrm{eff}}`` in the action
functionals.
"""
function noise_strength(sys::CoupledSDEs)
    Σ = covariance_matrix(sys)
    return sqrt(tr(Σ) / size(Σ, 1))
end
