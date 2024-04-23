"""
    fw_action(sys::CoupledSDEs, path, time; kwargs...)
Calculates the Freidlin-Wentzell action of a given `path` with time points `time` in a
drift field specified by the deterministic dynamics of `sys`.

The path must be a `(D x N)` matrix, where `D` is the dimensionality of the system `sys` and
`N` is the number of path points. The `time` array must have length `N`.

Returns a single number, which is the value of the action functional

``S_T[\\phi_t] = \\frac{1}{2} \\int_0^T || \\dot \\phi_t - b(\\phi_t) ||^2_Q dt``

where ``\\phi_t`` denotes the path in state space, ``b`` the drift field, and ``T`` the
total time of the path. The subscript ``Q`` refers to the
generalized norm ``||a||_Q^2 := \\langle a, Q^{-1} b \\rangle`` (see [`anorm`](@ref)). Here
``Q`` is the noise covariance matrix `sys.Σ`.

## Keyword arguments
* `cov_inv = nothing`: Inverse of the covariance matrix ``\\Sigma``.
  If `nothing`, ``\\Sigma^{-1}`` is computed automatically.
"""
function fw_action(sys::CoupledSDEs, path, time; cov_inv=nothing)
    all(diff(time) .≈ diff(time[1:2])) && error("Fw_action is only defined for equispaced time")
    # Inverse of covariance matrix
    A = isnothing(cov_inv) ? inv(covariance_matrix(sys)) : cov_inv

    # Compute action integral
    integrand = fw_integrand(sys, path, time, A)

    S = 0
    for i in 1:(size(path, 2) - 1)
        S += (integrand[i+1] + integrand[i])/2 * (time[i+1]-time[i])
    end
    S/2
end;

"""
    om_action(sys::CoupledSDEs, path, time; kwargs...)
Calculates the Onsager-Machlup action of a given `path` with time points `time` for the
drift field `sys.f` at noise strength `sys.σ`.

The path must be a `(D x N)` matrix, where `D` is the dimensionality of the system `sys` and
`N` is the number of path points. The `time` array must have length `N`.

Returns a single number, which is the value of the action functional

``I^{\\sigma}_T[\\phi_t] = \\frac{1}{2} \\int_0^T \\left( || \\dot \\phi_t - b(\\phi_t) ||^2_Q +
\\frac{\\sigma^2}{2} \\div(b) \\right) \\, dt``

where ``\\phi_t`` denotes the path in state space, ``b`` the drift field, ``T`` the total
time of the path, and ``\\sigma`` the noise strength. The subscript ``Q`` refers to the
generalized norm ``||a||_Q^2 := \\langle a, Q^{-1} b \\rangle`` (see [`anorm`](@ref)). Here
``Q`` is the noise covariance matrix `sys.Σ`.

## Keyword arguments
* `cov_inv = nothing`: Inverse of the covariance matrix ``\\Sigma``.
  If `nothing`, ``\\Sigma^{-1}`` is computed automatically.
"""
# function om_action(sys::CoupledSDEs, path, time; cov_inv=nothing)
#     # Inverse of covariance matrix
#     (cov_inv == nothing) ? A = inv(sys.Σ) : A = cov_inv

#     # Compute action integral
#     S = 0
#     for i in 1:(size(path, 2) - 1)
#         S += sys.σ^2/2 * ((div_b(sys, path[:,i+1]) + div_b(sys, path[:,i]))/2 *
#             (time[i+1]-time[i]))
#     end
#     fw_action(sys, path, time, cov_inv=A) + S/2
# end;

"""
    action(sys::CoupledSDEs, path::Matrix, time, functional; kwargs...)
Computes the action functional specified by `functional` for a given CoupledSDEs `sys` and
`path` parameterized by `time`.

* `functional = "FW"`: Returns the Freidlin-Wentzell action ([`fw_action`](@ref))
* `functional = "OM"`: Returns the Onsager-Machlup action ([`om_action`](@ref))
"""
function action(sys::CoupledSDEs, path::Matrix, time, functional; kwargs...)
    all(diff(time) .≈ diff(time[1:2])) && error("action is only defined for equispaced time")
    if functional == "FW"
        return fw_action(sys, path, time; kwargs...)
    elseif functional == "OM"
        return om_action(sys, path, time; kwargs...)
    end
end;

"""
    geometric_action(sys::CoupledSDEs, path, arclength=1; kwargs...)
Calculates the geometric action of a given `path` with specified `arclength` for the drift
field `sys.f`.

For a given path ``\\varphi``, the geometric action ``\\bar S`` corresponds to the minimum
of the Freidlin-Wentzell action ``S_T(\\phi)`` over all travel times ``T>0``, where ``\\phi``
denotes the path's parameterization in physical time (see [`fw_action`](@ref)). It is given
by the integral

``\\bar S[\\varphi] = \\int_0^L \\left( ||\\varphi'||_Q \\, ||b(\\varphi)||_Q - \\langle \\varphi', \\,
    b(\\varphi) \\rangle_Q \\right) \\, ds``

where ``s`` is the arclength coordinate, ``L`` the arclength, ``b`` the drift field, and the
subscript ``Q`` refers to the generalized dot product ``\\langle a, b \\rangle_Q := a^{\\top}
\\cdot Q^{-1} b`` (see [`anorm`](@ref)). Here ``Q`` is the noise covariance matrix `sys.Σ`.

## Keyword arguments
* `cov_inv = nothing`: Inverse of the covariance matrix ``\\Sigma``.
  If `nothing`, ``\\Sigma^{-1}`` is computed automatically.

Returns the value of the geometric action ``\\bar S``.
"""
function geometric_action(sys::CoupledSDEs, path, arclength=1.0; cov_inv=nothing)
    N = size(path, 2)
    v = path_velocity(path, range(0, arclength, length=N), order=4)
    A = isnothing(cov_inv) ? inv(covariance_matrix(sys)) : cov_inv

    b(x) = drift(sys, x)

    integrand = zeros(N)
    for i in 1:N
        drift = b(path[:,i])
        integrand[i] = anorm(v[:,i], A)*anorm(drift, A) - dot(v[:,i], A, drift)
    end

    S = 0
    for i in 1:N-1
        S += (integrand[i+1] + integrand[i])/2
    end
    S * arclength/(N-1)
end

"""
    fw_integrand(sys::CoupledSDEs, path, time, A)
Computes the squared ``A``-norm ``|| \\dot \\phi_t - b ||^2_A`` (see `fw_action` for
details). Returns a vector of length `N` containing the values of the above squared norm for
each time point in the vector `time`.
"""
function fw_integrand(sys::CoupledSDEs, path, time, A)
    v = path_velocity(path, time, order=4)
    sqnorm = zeros(size(path, 2))
    b(x) = drift(sys, x)
    for i in 1:size(path, 2)
        # assumes the drift is time independent
        drift = b(path[:,i])
        sqnorm[i] = anorm(v[:,i] - drift, A, square=true)
    end
    sqnorm
end;

"""
    div_b(sys::CoupledSDEs, x)
Computes the divergence of the drift field `sys.f` at the given point `x`.
"""
function div_b(sys::CoupledSDEs, x)
    b(x) = drift(sys, x)
    tr(ForwardDiff.jacobian(b, x))
end;

"""
    path_velocity(path, time; order=4)
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
        v[:,1] .= (path[:,2] .- path[:,1])/(time[2] - time[1])
        v[:,end] .= (path[:,end] .- path[:,end-1])/(time[end] - time[end-1])
        # 2nd order central differences for internal points
        for i in 2:(size(path, 2) - 1)
            v[:,i] .= (path[:,i+1] .- path[:,i-1])/(time[i+1] - time[i-1])
        end

    elseif order == 4
        # 1st order forward/backward differences for end points
        v[:,1] .= (path[:,2] .- path[:,1])/(time[2] - time[1])
        v[:,end] .= (path[:,end] .- path[:,end-1])/(time[end] - time[end-1])
        # 2nd order central differences for neighbors of end points
        v[:,2] .= (path[:,3] .- path[:,1])/(time[3] - time[1])
        v[:,end-1] .= (path[:,end] .- path[:,end-2])/(time[end] - time[end-2])
        # 4th order central differences for internal points
        for i in 3:(size(path, 2) - 2)
            v[:,i] .= ((-path[:,i+2] .+ 8*path[:,i+1] .- 8*path[:,i-1] .+ path[:,i-2])/
                (6*(time[i+1] - time[i-1])))
        end
    end
    v
end;
