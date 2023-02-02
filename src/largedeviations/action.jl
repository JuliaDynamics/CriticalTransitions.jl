include("../StochSystem.jl")
include("../utils.jl")

"""
    fw_action(sys::StochSystem, path, time)
Calculates the Freidlin-Wentzell action of a given `path` with time points `time` in a 
drift field specified by the deterministic dynamics of `sys`.

The path must be a `(D x N)` matrix, where `D` is the dimensionality of the system `sys` and
`N` is the number of path points. The `time` array must have length `N`.

Returns a single number, which is the value of the action functional

``S_T[\\phi_t] = \\frac{1}{2} \\int_0^T || \\dot \\phi_t - b ||^2_A dt``

where ``\\phi_t`` denotes the path in state space, ``b`` the drift field, and ``T`` the
total time of the path. The norm is taken with respect to the matrix ``A``
(see [`anorm`](@ref)), which is the Cholesky decomposition of the covariance matrix
``\\Sigma = AA^\\top`` specified by `sys.Σ`.
"""
function fw_action(sys::StochSystem, path, time)
    integrand = fw_integrand(sys, path, time)
    S = 0
    for i in 1:(size(path, 2) - 1)
        S += (integrand[i+1] + integrand[i])/2 * (time[i+1]-time[i])
    end
    S/2
end;

"""
    om_action(sys::StochSystem, path, time)
Calculates the Onsager-Machlup action of a given `path` with time points `time` for the
drift field `sys.f` at noise strength `sys.σ`.

The path must be a `(D x N)` matrix, where `D` is the dimensionality of the system `sys` and
`N` is the number of path points. The `time` array must have length `N`.

Returns a single number, which is the value of the action functional

``I^{\\sigma}_T[\\phi_t] = \\frac{1}{2} \\int_0^T || \\dot \\phi_t - b ||^2_A +
\\frac{\\sigma^2}{2} \\nabla \\dot b dt``

where ``\\phi_t`` denotes the path in state space, ``b`` the drift field, ``T`` the total
time of the path, and ``\\sigma`` the noise strength. The norm is taken with respect to the
matrix ``A`` (see [anorm](@ref)), which is the Cholesky decomposition of the covariance
matrix ``\\Sigma = AA^\\top`` specified by `sys.Σ`.
"""
function om_action(sys::StochSystem, path, time)
    S = 0
    for i in 1:(size(path, 2) - 1)
        S += sys.σ^2/2 * ((div_b(sys, path[:,i+1]) + div_b(sys, path[:,i]))/2 *
            (time[i+1]-time[i]))
    end
    fw_action(sys, path, time) + S/2
end;

"""
    fw_integrand(sys::StochSystem, path, time)
Computes the squared ``A``-norm ``|| \\dot \\phi_t - b ||^2_A`` (see `fw_action` for
details). Returns a vector of length `N` containing the values of the above squared norm for
each time point in the vector `time`.
"""
function fw_integrand(sys::StochSystem, path, time)
    v = path_velocity(path, time, order=4)
    A = inv(sys.Σ)
    sqnorm = zeros(size(path, 2))
    for i in 1:size(path, 2)
        drift = sys.f(path[:,i], p(sys), time[i])
        sqnorm[i] = anorm(v[:,i] - drift, A, square=true)
    end
    sqnorm
end;

"""
    div_b(sys::StochSystem, x)
Computes the divergence of the drift field `sys.f` at the given point `x`.
"""
function div_b(sys::StochSystem, x)
    b(x) = sys.f(x, p(sys), 0)
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