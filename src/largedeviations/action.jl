"""
$(TYPEDSIGNATURES)

Calculates the Freidlin-Wentzell action of a given `path` with time points `time` in a
drift field specified by the deterministic dynamics `f = dynamic_rule(sys)` and
(normalized) noise covariance matrix `covariance_matrix(sys)`.

The path must be a `(D x N)` matrix, where `D` is the dimensionality of the system `sys` and
`N` is the number of path points. The `time` array must have length `N`.

Returns a single number, which is the value of the action functional

``S_T[\\phi_t] = \\frac{1}{2} \\int_0^T || \\dot \\phi_t - f(\\phi_t) ||^2_Q \\text{d}t``

where ``\\phi_t`` denotes the path in state space, ``b`` the drift field, and ``T`` the
total time of the path. The subscript ``Q`` refers to the
generalized norm ``||a||_Q^2 := \\langle a, Q^{-1} b \\rangle`` (see `anorm`). Here
``Q`` is the noise covariance matrix normalized by ``D/L_1(Q)``, with ``L_1`` being the
L1 matrix norm.
"""
function fw_action(sys::CoupledSDEs, path, time)
    @assert all(diff(time) .≈ diff(time[1:2])) "Freidlin-Wentzell action is only defined for equispaced time"
    # Inverse of covariance matrix
    A = inv(normalize_covariance!(covariance_matrix(sys)))

    # Compute action integral
    integrand = fw_integrand(sys, path, time, A)

    S = 0
    for i in 1:(size(path, 2) - 1)
        S += (integrand[i + 1] + integrand[i]) / 2 * (time[i + 1] - time[i])
    end
    return S / 2
end;

"""
    om_action(sys::CoupledSDEs, path, time, noise_strength)

Calculates the Onsager-Machlup action of a given `path` with time points `time` for the drift field `f = dynamic_rule(sys)` at given `noise_strength`.

The path must be a `(D x N)` matrix, where `D` is the dimensionality of the system `sys` and
`N` is the number of path points. The `time` array must have length `N`.

Returns a single number, which is the value of the action functional

``S^{\\sigma}_T[\\phi_t] = \\frac{1}{2} \\int_0^T \\left( || \\dot \\phi_t - f(\\phi_t) ||^2_Q +
\\sigma^2 \\nabla \\cdot f \\right) \\, \\text{d} t``

where ``\\phi_t`` denotes the path in state space, ``b`` the drift field, ``T`` the total
time of the path, and ``\\sigma`` the noise strength. The subscript ``Q`` refers to the
generalized norm ``||a||_Q^2 := \\langle a, Q^{-1} b \\rangle`` (see `anorm`). Here
``Q`` is the noise covariance matrix normalized by ``D/L_1(Q)``, with ``L_1`` being the
L1 matrix norm.
"""
function om_action(sys::CoupledSDEs, path, time, noise_strength)
    @assert all(diff(time) .≈ diff(time[1:2])) "Onsager-Machlup action is only defined for equispaced time"

    σ = noise_strength
    # Compute action integral
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
function action(sys::CoupledSDEs, path::Matrix, time, functional; noise_strength=nothing)
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

Calculates the geometric action of a given `path` with specified `arclength` for the drift
field specified by the deterministic dynamics `f = dynamic_rule(sys)` and
(normalized) noise covariance matrix `covariance_matrix(sys)`.

For a given path ``\\varphi``, the geometric action ``\\bar S`` corresponds to the minimum
of the Freidlin-Wentzell action ``S_T(\\varphi)`` over all travel times ``T>0``, where ``\\varphi``
denotes the path's parameterization in physical time (see [`fw_action`](@ref)). It is given
by the integral

``\\bar S[\\varphi] = \\int_0^L \\left( ||\\varphi'||_Q \\, ||f(\\varphi)||_Q - \\langle \\varphi', \\,
    f(\\varphi) \\rangle_Q \\right) \\, \\text{d}s``

where ``s`` is the arclength coordinate, ``L`` the arclength, ``f`` the drift field, and the
subscript ``Q`` refers to the generalized dot product ``\\langle a, b \\rangle_Q := a^{\\top}
\\cdot Q^{-1} b`` (see `anorm`). Here
``Q`` is the noise covariance matrix normalized by ``D/L_1(Q)``, with ``L_1`` being the
L1 matrix norm.

Returns the value of the geometric action ``\\bar S``.
"""
function geometric_action(sys::CoupledSDEs, path, arclength=1.0)
    A = inv(normalize_covariance!(covariance_matrix(sys)))
    b(x) = drift(sys, x)
    return _geometric_action_from_drift(b, path, arclength, A)
end

"""
    geometric_action(b::Function, path, arclength=1.0; A=nothing)

Geometric Freidlin-Wentzell action for a drift field `b` and a discrete `path`.

This is a drift-only convenience overload that uses the same integrand as
`geometric_action(sys::CoupledSDEs, ...)`, but requires an explicit inverse covariance
(`A = Q^{-1}`) if you want anything other than the identity metric.

If `A === nothing`, it defaults to the identity metric.

The `path` must be a `(D × N)` matrix with points in columns.
"""
function geometric_action(b::Function, path, arclength=1.0; A=nothing)
    if A === nothing
        D = size(path, 1)
        A = LinearAlgebra.Diagonal(ones(eltype(path), D))
    end
    return _geometric_action_from_drift(b, path, arclength, A)
end

function _geometric_action_from_drift(b::Function, path, arclength::Real, A)
    N = size(path, 2)
    v = path_velocity(path, range(0, arclength; length=N); order=4)

    integrand = zeros(eltype(path), N)
    for i in 1:N
        drift = b(path[:, i])
        integrand[i] = anorm(v[:, i], A) * anorm(drift, A) - dot(v[:, i], A, drift)
    end

    S = zero(eltype(path))
    for i in 1:(N - 1)
        S += (integrand[i + 1] + integrand[i]) / 2
    end
    return S * arclength / (N - 1)
end

"""
$(TYPEDSIGNATURES)

Computes the squared ``A``-norm ``|| \\dot \\phi_t - b ||^2_A`` (see `fw_action` for
details). Returns a vector of length `N` containing the values of the above squared norm for each time point in the vector `time`.
"""
function fw_integrand(sys::CoupledSDEs, path, time, A)
    v = path_velocity(path, time; order=4)
    sqnorm = zeros(size(path, 2))
    b(x) = drift(sys, x)
    for i in axes(path, 2)
        # assumes the drift is time independent
        drift = b(path[:, i])
        sqnorm[i] = anorm(v[:, i] - drift, A; square=true)
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
function path_velocity(path, time; order=4)
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
