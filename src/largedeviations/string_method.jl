"""
    $(TYPEDSIGNATURES)

Compute the string method for a given system using [E et al. (2007)](https://doi.org/10.1063/1.2720838).

The string method is an iterative algorithm used to find minimum energy path (MEP) between two points in phase space. It works by evolving a discretized path (string) according to the system's drift while maintaining equal arc-length parametrization between points.

This implementation allows for computation between arbitrary points, not just stable fixed points.

# Arguments
- `sys::SgmamSystem`: The doubled phase space system for which the string method is computed
- `x_initial`: Initial path discretized as a matrix where each column represents a point on the path
- `ϵ::Float64`: Step size for the evolution step
- `iterations::Int64`: Maximum number of iterations for path convergence
- `show_progress::Bool`: Whether to display a progress meter during computation

# Returns
- `x`: The final converged path representing the MEP
"""
function string_method(
    sys::Union{SgmamSystem,Function},
    x_initial::Matrix;
    ϵ::Float64=1e-1,
    iterations::Int64=1000,
    show_progress::Bool=false,
)
    Nx, Nt = size(x_initial)
    s = range(0; stop=1, length=Nt)
    x, alpha = init_allocation_string(x_initial, Nt)

    progress = Progress(iterations; dt=0.5, enabled=show_progress)
    for i in 1:iterations
        update_x!(x, sys, ϵ)
        # reset initial and final points to allow for string computation
        # between points that are not stable fixed points
        x[:, 1] = x_initial[:, 1]
        x[:, end] = x_initial[:, end]

        # reparametrize to arclength
        interpolate_path!(x, alpha, s)

        next!(progress; showvalues=[("iterations", i)])
    end
    return StateSpaceSet(x')
end
"""
    $(TYPEDSIGNATURES)

Compute the string method for a given system using [E et al. (2007)](https://doi.org/10.1063/1.2720838).

The string method is an iterative algorithm used to find minimum energy path (MEP) between two points in phase space. It works by evolving a discretized path (string) according to the system's drift while maintaining equal arc-length parametrization between points.

This implementation allows for computation between arbitrary points, not just stable fixed points.

# Arguments
- `sys::CoupledSDEs`: The system for which the string method is computed
- `x_initial`: Initial path discretized as a matrix where each column represents a point on the path
- `ϵ::Float64`: Step size for the evolution step
- `iterations::Int64`: Maximum number of iterations for path convergence
- `show_progress::Bool`: Whether to display a progress meter during computation

# Returns
- `x`: The final converged path representing the MEP
"""
function string_method(sys::CoupledSDEs, init; kwargs...)
    b(x) = drift(sys, x)
    return string_method(b, init; kwargs...)
end

function string_method(
    b::Union{SgmamSystem,Function},
    x_initial::StateSpaceSet{D};
    ϵ::Float64=1e-1,
    iterations::Int64=1000,
    show_progress::Bool=false,
) where {D}
    Nt = length(x_initial)
    s = range(0; stop=1, length=Nt)
    x, alpha = init_allocation_string(x_initial, Nt)

    progress = Progress(iterations; dt=0.5, enabled=show_progress)
    for i in 1:iterations
        update_x!(x, b, ϵ)

        # reparametrize to arclength
        x = interpolate_path(x, alpha, s)

        next!(progress; showvalues=[("iterations", i)])
    end
    return x
end

function update_x!(x::Matrix, sys::SgmamSystem, ϵ::Float64)
    return x += ϵ * sys.H_p(x, 0 * x) # euler integration
end
function update_x!(x::StateSpaceSet, sys::SgmamSystem, ϵ::Float64)
    return x += ϵ .* vec(sys.H_p(x, 0 * Matrix(x))) # euler integration
end
function update_x!(x::StateSpaceSet, b::Function, ϵ::Float64)
    # do not touch the initial and final points to allow for string computation
    # between points that are not stable fixed points
    for j in 2:(length(x) - 1)
        x[j] += ϵ * b(x[j]) # euler integration
    end
end
function update_x!(x::Matrix, b::Function, ϵ::Float64)
    # do not touch the initial and final points to allow for string computation
    # between points that are not stable fixed points
    for j in 2:(size(x)[2] - 1)
        x[:,j] += ϵ * b(x[:,j]) # euler integration
    end
end

function init_allocation_string(x_initial, Nt)
    # preallocate
    x = deepcopy(x_initial)
    alpha = zeros(Nt)
    return x, alpha
end
