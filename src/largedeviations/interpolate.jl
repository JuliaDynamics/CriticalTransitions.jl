"""
    interpolate_path!(path, α, s)

Interpolate a path in-place to ensure uniform spacing between points.

# Arguments
- `path`: Matrix of size (D, N) containing the path points, where D is the dimension and N is the number of points
- `α`: Vector of length N to store the normalized cumulative distances
- `s`: Vector of length N containing the desired interpolation points (typically uniform from 0 to 1)

# Details
The function performs these steps:
1. Computes distances between consecutive points
2. Normalizes cumulative distances to [0,1] interval
3. Interpolates each dimension of the path using the normalized distances

The interpolation is performed in-place, modifying both `path` and `α`.
"""
function interpolate_path!(path::Matrix, α, s)
    α[2:end] .= vec(sqrt.(sum(diff(path; dims=2) .^ 2; dims=1)))
    α .= cumsum(α; dims=1)
    α .= α ./ last(α)
    for dof in 1:size(path, 1)
        path[dof, :] .= linear_interpolation(α, path[dof, :])(s)
    end
    return nothing
end
function interpolate_path(path::StateSpaceSet{D}, α, s) where {D}
    Matrix
    α[2:end] .= vec(sqrt.(sum.(map(x -> x .^ 2, diff(path)))))
    α .= cumsum(α; dims=1)
    α .= α ./ last(α)
    return StateSpaceSet([linear_interpolation(α, path[:, dof])(s) for dof in 1:D]...)
end

"""
$(TYPEDSIGNATURES)

Perform an in-place Euler step

### Fields
* `u` - Initial path
* `sys` - doubled phase space system
* `ϵ` - Time step
"""
function updateEuler!(x::Matrix, sys::SgmamSystem, ϵ::Float64)
    return x += ϵ * sys.H_p(x, 0 * x) # euler integration
end
"""
$(TYPEDSIGNATURES)

Perform an in-place Euler step

### Fields
* `u` - Initial path
* `sys` - doubled phase space system
* `ϵ` - Time step
"""
function updateEuler!(x::StateSpaceSet, sys::SgmamSystem, ϵ::Float64)
    return x += ϵ .* vec(sys.H_p(x, 0 * Matrix(x))) # euler integration
end
"""
$(TYPEDSIGNATURES)

Perform an in-place Euler step

### Fields
* `u` - Initial path
* `b` - Gradient of energy
* `ϵ` - Time step
"""
function updateEuler!(x::StateSpaceSet, b::Function, ϵ::Float64)
    # do not touch the initial and final points to allow for string computation
    # between points that are not stable fixed points
    for j in 2:(length(x) - 1)
        x[j] += ϵ * b(x[j]) # euler integration
    end
end
"""
$(TYPEDSIGNATURES)

Perform an in-place Euler step

### Fields
* `u` - Initial path
* `b` - Gradient of energy
* `ϵ` - Time step
"""
function updateEuler!(x::Matrix, b::Function, ϵ::Float64)
    # do not touch the initial and final points to allow for string computation
    # between points that are not stable fixed points
    for j in 2:(size(x)[2] - 1)
        x[:, j] += ϵ * b(x[:, j]) # euler integration
    end
end

"""
`stepRK4!(u, b, ϵ)`: Perform an in place RK4 step

### Fields
* `u` - Initial state
* `b` - Gradient of energy
* `ϵ` - Time step
"""
function stepRK4!(u, b, ϵ)
    gradV1 = b(u)
    gradV2 = b(u - 0.5 * ϵ * gradV1)
    gradV3 = b(u - 0.5 * ϵ * gradV2)
    gradV4 = b(u - ϵ * gradV3)

    @. u = u - ϵ / 6 * (gradV1 + 2 * gradV2 + 2 * gradV3 + gradV4)

    return u
end
