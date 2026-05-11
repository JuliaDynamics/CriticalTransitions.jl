"""
$(TYPEDSIGNATURES)

Calculates the Freidlin-Wentzell action of a given `path` with time points `time` in
the drift field `f = dynamic_rule(sys)` and (normalized) diffusion tensor of `sys`.

The path must be a `(D × N)` matrix where `D` is the dimensionality of `sys` and `N`
is the number of path points. `time` must have length `N`.

Returns a single number, the value of the action functional

``S_T[\\phi_t] = \\frac{1}{2} \\int_0^T \\lVert \\dot \\phi_t - f(\\phi_t) \\rVert
    ^2_{a(\\phi_t)^{-1}} \\, \\mathrm{d}t``

where ``a(x) = \\sigma(x)\\sigma(x)^\\top`` is the diffusion tensor and the subscript
``a^{-1}`` denotes the generalized norm
``\\lVert v\\rVert^2_{a^{-1}} := \\langle v,\\, a^{-1} v\\rangle`` (see [`anorm`](@ref)).

This is the Freidlin-Wentzell **rate function**: it is invariant under uniform rescaling
of the noise. The diffusion is normalized by a single scale extracted from `sys`,
matching the long-standing CT convention; see the Large Deviations manual page for
details.
"""
function fw_action(sys::CoupledSDEs, path, time)
    @assert all(diff(time) .≈ diff(time[1:2])) "Freidlin-Wentzell action is only defined for equispaced time"
    proper_MAM_system(sys)
    A_at = _inv_covariance_at(sys)
    integrand = fw_integrand(sys, path, time, A_at)

    S = 0
    for i in 1:(size(path, 2) - 1)
        S += (integrand[i + 1] + integrand[i]) / 2 * (time[i + 1] - time[i])
    end
    return S / 2
end;

"""
    om_action(sys::CoupledSDEs, path, time, noise_strength)

Calculates the Onsager-Machlup action of a given `path` with time points `time` for the
drift `f = dynamic_rule(sys)` at given `noise_strength` ``\\sigma``.

Returns a single number,

``S^{\\sigma}_T[\\phi_t] = \\frac{1}{2} \\int_0^T \\left( \\lVert \\dot \\phi_t -
    f(\\phi_t) \\rVert^2_{a^{-1}} + \\sigma^2 \\nabla \\cdot f \\right) \\, \\mathrm{d}t``

The ``\\sigma^2 \\nabla\\cdot f`` divergence correction assumes ``a`` is **constant**,
i.e. additive noise. Multiplicative noise is rejected with an `ArgumentError`; for
state-dependent ``a(x)`` use [`fw_action`](@ref) (which has no divergence correction).
"""
function om_action(sys::CoupledSDEs, path, time, noise_strength)
    @assert all(diff(time) .≈ diff(time[1:2])) "Onsager-Machlup action is only defined for equispaced time"
    sys.noise_type[:additive] || throw(
        ArgumentError(
            "om_action is only defined for additive noise: the σ²∇·f divergence " *
                "correction is invalid when a(x) is state-dependent. Use fw_action instead.",
        ),
    )

    σ = noise_strength
    S = 0.0
    for i in 1:(size(path, 2) - 1)
        S +=
            σ^2 / 2 * (
            (div_drift(sys, path[:, i + 1]) + div_drift(sys, path[:, i])) / 2 *
                (time[i + 1] - time[i])
        )
    end
    return fw_action(sys, path, time) + S / 2
end;

"""
$(TYPEDSIGNATURES)

Computes the action functional specified by `functional` for a given CoupledSDEs `sys` and
`path` parameterized by `time`.

* `functional = "FW"`: Returns the Freidlin-Wentzell action ([`fw_action`](@ref))
* `functional = "OM"`: Returns the Onsager-Machlup action ([`om_action`](@ref))
"""
function action(sys::CoupledSDEs, path::Matrix, time, functional; noise_strength = nothing)
    S = 0.0
    if functional == "FW"
        S = fw_action(sys, path, time)
    elseif functional == "OM"
        S = om_action(sys, path, time, noise_strength)
    end
    return S
end;

"""
$(TYPEDSIGNATURES)

Calculates the geometric Freidlin-Wentzell action of a given `path` with specified
`arclength` for the drift `f = dynamic_rule(sys)` and (normalized) diffusion tensor of
`sys`.

For a given path ``\\varphi``, the geometric action ``\\bar S`` corresponds to the
minimum of the Freidlin-Wentzell action ``S_T(\\varphi)`` over all travel times
``T > 0`` (see [`fw_action`](@ref)). It is given by

``\\bar S[\\varphi] = \\int_0^L \\left( \\lVert\\varphi'\\rVert_{a^{-1}}\\,
    \\lVert f(\\varphi)\\rVert_{a^{-1}} - \\langle \\varphi',\\, f(\\varphi)\\rangle_{a^{-1}}
    \\right) \\, \\mathrm{d}s``

where ``s`` is the arclength coordinate, ``L`` the arclength, and the subscript
``a^{-1}`` denotes the generalized norm/inner product (see [`anorm`](@ref)). The same
noise normalization convention as [`fw_action`](@ref) applies.
"""
function geometric_action(sys::CoupledSDEs, path, arclength = 1.0)
    proper_MAM_system(sys)
    A_at = _inv_covariance_at(sys)
    b(x) = drift(sys, x)
    return _geometric_action_from_drift(b, path, arclength, A_at)
end

"""
    geometric_action(b::Function, path, arclength=1.0; A=nothing)

Geometric Freidlin-Wentzell action for a drift field `b` and a discrete `path`. The
metric argument `A` selects how the inverse covariance is evaluated at each point:

* `A = nothing` (default): identity inverse covariance everywhere.
* `A::AbstractMatrix`: constant inverse covariance ``a^{-1}``, used at every path point.
* `A::Function`: callable returning ``a(x)`` (the uninverted covariance) at state `x`;
  inverted internally at each path point.

The `path` must be a `(D × N)` matrix with points in columns.
"""
function geometric_action(b::Function, path, arclength = 1.0; A = nothing)
    A_at = if A === nothing
        D = size(path, 1)
        Returns(LinearAlgebra.Diagonal(ones(eltype(path), D)))
    elseif A isa AbstractMatrix
        Returns(A)
    else
        let A = A
            x -> inv(A(x))
        end
    end
    return _geometric_action_from_drift(b, path, arclength, A_at)
end

# Build a callable `x -> a_shape(x)^{-1}` for the given system, where
# `a_shape(x) = a(x) / s` and `s = L_1(a(x_0))/D` is the noise-magnitude scale extracted
# at the reference state `x_0 = current_state(sys)`. See the Large Deviations manual page
# ("Action normalisation") for the rationale.
function _inv_covariance_at(sys::CoupledSDEs)
    if sys.noise_type[:additive]
        a = covariance_matrix(sys)
        s = LinearAlgebra.norm(a, 1) / size(a, 1)
        return Returns(s * inv(a))  # = inv(a/s) = a_shape⁻¹
    end
    g = diffusion_function(sys)
    params = current_parameters(sys)
    σ_ref = g(current_state(sys), params, 0.0)
    a_ref = σ_ref * σ_ref'
    s = LinearAlgebra.norm(a_ref, 1) / size(a_ref, 1)
    return let g = g, params = params, s = s
        x -> begin
            σ = g(x, params, 0.0)
            s * inv(σ * σ')
        end
    end
end

function _geometric_action_from_drift(b::Function, path, arclength::Real, A_at)
    N = size(path, 2)
    v = path_velocity(path, range(0, arclength; length = N); order = 4)

    integrand = zeros(eltype(path), N)
    for i in 1:N
        xi = view(path, :, i)
        Ai = A_at(xi)
        drift_i = b(xi)
        integrand[i] = anorm(v[:, i], Ai) * anorm(drift_i, Ai) - dot(v[:, i], Ai, drift_i)
    end

    S = zero(eltype(path))
    for i in 1:(N - 1)
        S += (integrand[i + 1] + integrand[i]) / 2
    end
    return S * arclength / (N - 1)
end

"""
$(TYPEDSIGNATURES)

Computes the squared ``a^{-1}``-norm ``\\lVert \\dot \\phi_t - b\\rVert^2_{a(\\phi_t)^{-1}}``
at each path point (see [`fw_action`](@ref) for details). `A` may be either a constant
`AbstractMatrix` (used as the inverse diffusion tensor at every path point) or a callable
`x -> a(x)^{-1}` for state-dependent diffusion. Returns a vector of length `N`.
"""
function fw_integrand(sys::CoupledSDEs, path, time, A)
    v = path_velocity(path, time; order = 4)
    sqnorm = zeros(size(path, 2))
    for i in axes(path, 2)
        xi = view(path, :, i)
        drift_i = drift(sys, xi)
        Ai = A isa AbstractMatrix ? A : A(xi)
        sqnorm[i] = anorm(v[:, i] - drift_i, Ai; square = true)
    end
    return sqnorm
end;

"""
$(TYPEDSIGNATURES)

Returns the velocity along a given `path` with time points given by `time`.

## Keyword arguments
* `order = 4`: Accuracy of the finite difference approximation.
  `4`th order corresponds to a five-point stencil, `2`nd order to a three-point stencil.
  In both cases, central differences are used except at the end points, where a forward or
  backward difference is used.
"""
function path_velocity(path, time; order = 4)
    v = zeros(size(path))

    if order == 2
        # 1st order forward/backward differences for end points
        v[:, 1] .= (path[:, 2] .- path[:, 1]) / (time[2] - time[1])
        v[:, end] .= (path[:, end] .- path[:, end - 1]) / (time[end] - time[end - 1])
        # 2nd order central differences for internal points
        for i in 2:(size(path, 2) - 1)
            v[:, i] .= (path[:, i + 1] .- path[:, i - 1]) / (time[i + 1] - time[i - 1])
        end

    elseif order == 4
        # 1st order forward/backward differences for end points
        v[:, 1] .= (path[:, 2] .- path[:, 1]) / (time[2] - time[1])
        v[:, end] .= (path[:, end] .- path[:, end - 1]) / (time[end] - time[end - 1])
        # 2nd order central differences for neighbors of end points
        v[:, 2] .= (path[:, 3] .- path[:, 1]) / (time[3] - time[1])
        v[:, end - 1] .= (path[:, end] .- path[:, end - 2]) / (time[end] - time[end - 2])
        # 4th order central differences for internal points
        for i in 3:(size(path, 2) - 2)
            v[:, i] .= (
                (
                    -path[:, i + 2] .+ 8 * path[:, i + 1] .- 8 * path[:, i - 1] .+
                        path[:, i - 2]
                ) / (6 * (time[i + 1] - time[i - 1]))
            )
        end
    end
    return v
end;
