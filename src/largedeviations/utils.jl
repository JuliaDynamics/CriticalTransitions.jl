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
function interpolate_path(path::StateSpaceSet{D}, α, s) where D
    Matrix
    α[2:end] .= vec(sqrt.(sum.(map(x -> x .^ 2, diff(path)))))
    α .= cumsum(α; dims=1)
    α .= α ./ last(α)
    return StateSpaceSet([linear_interpolation(α, path[:,dof])(s) for dof in 1:D]...)
end

"""
`stepEuler!(u, b, Δt)`: Perform an in place Euler step

### Fields
* `u` - Initial state
* `b` - Gradient of energy
* `Δt` - Time step
"""
function stepEuler!(u, b, Δt)
    gradV = b(u)
    @. u = u - Δt * gradV

    return u
end

"""
`stepRK4!(u, b, Δt)`: Perform an in place RK4 step

### Fields
* `u` - Initial state
* `b` - Gradient of energy
* `Δt` - Time step
"""
function stepRK4!(u, b, Δt)
    gradV1 = b(u)
    gradV2 = b(u - 0.5 * Δt * gradV1)
    gradV3 = b(u - 0.5 * Δt * gradV2)
    gradV4 = b(u - Δt * gradV3)

    @. u = u - Δt / 6 * (gradV1 + 2 * gradV2 + 2 * gradV3 + gradV4)

    return u
end
