include("../StochSystem.jl")
include("../utils.jl")

"""
    fw_action(sys::StochSystem, path, time)
Calculates the Freidlin-Wentzell action of a given `path` with time points `time` in a 
drift field specified by the deterministic dynamics of `sys`.

The path must be a `(D x N)` matrix, where `D` is the dimensionality of the system `sys` and
`N` is the number of path points. The `time` array must have length `N`.

Returns a single number, which is the value of the action integral

``S_T[\\phi_t] = \\frac{1}{2} \\int_0^T || \\dot \\phi_t - b ||^2 dt``

where ``\\phi_t`` denotes the path and ``b`` the drift field.
"""
function fw_action(sys::StochSystem, path, time)
    integrand = fw_integrand(sys, path, time)
    S = 0
    for i in 1:(size(path, 2) - 1)
        S += (integrand[i+1] + integrand[i])/2 * (time[i+1]-time[i])
    end
    S/2
end;

function om_action(sys::StochSystem, path, time)
    integrand = fw_integrand(sys, path, time)
    S = 0
    for i in 1:(size(path, 2) - 1)
        S += (integrand[i+1] + integrand[i])/2 * (time[i+1]-time[i])
        S += sys.σ^2 * (div_b(sys, path[:,i+1]) + div_b(sys, path[:,i]))/2 * (time[i+1]-time[i])
    end
    S/2
end;

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

function div_b(sys::StochSystem, x)
    b(x) = sys.f(x, p(sys), 0)
    tr(ForwardDiff.jacobian(b, x))
end;

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