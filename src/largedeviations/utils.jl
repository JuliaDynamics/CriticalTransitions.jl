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
    α[2:end] .= vec(sqrt.(sum(diff(path; dims = 2) .^ 2; dims = 1)))
    α .= cumsum(α; dims = 1)
    α .= α ./ last(α)
    for dof in 1:size(path, 1)
        path[dof, :] .= linear_interpolation(α, path[dof, :])(s)
    end
    return nothing
end
function interpolate_path(path::StateSpaceSet{D}, α, s) where {D}
    Matrix
    α[2:end] .= vec(sqrt.(sum.(map(x -> x .^ 2, diff(path)))))
    α .= cumsum(α; dims = 1)
    α .= α ./ last(α)
    return StateSpaceSet([linear_interpolation(α, path[:, dof])(s) for dof in 1:D]...)
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

"""
    proper_MAM_system(ds::CoupledSDEs)

Validate that `ds` is admissible as input to the MAM / gMAM / sgMAM family.

All variants require `:autonomous` noise. For systems with additive noise, `:invertible`
must also hold (constant diffusion matrix must be invertible). For multiplicative /
state-dependent noise, invertibility of ``a(x)`` along the path is the caller's
responsibility.
"""
function proper_MAM_system(ds::CoupledSDEs)
    if !ds.noise_type[:autonomous]
        throw(
            ArgumentError(
                "The minimal action method is only applicable for autonomous noise; the system's noise is non-autonomous.",
            ),
        )
    end
    if ds.noise_type[:additive] && !ds.noise_type[:invertible]
        throw(
            ArgumentError(
                "The minimal action method requires an invertible constant diffusion matrix for additive noise.",
            ),
        )
    end
    return
end

# Helpers for the three-way dispatch on `a_func(x)` return types in the gMAM/sgMAM
# solvers and action integrands:
#
#   * `nothing`         → additive a = I (shared tridiagonal / shared metric)
#   * `AbstractVector`  → diagonal a(x) (per-DOF fast path)
#   * `Diagonal`        → diagonal a(x) wrapped as a LinearAlgebra.Diagonal; we accept it
#                         on the same fast path so hand-rolled `a_func(x) = Diagonal(...)`
#                         doesn't silently fall into the slow general-matrix code path.
#   * other `AbstractMatrix` → general state-dependent a(x) (sparse block-tridiagonal /
#                              dense inverse).
@inline _is_diagonal_a(a) = a isa AbstractVector || a isa LinearAlgebra.Diagonal
@inline _diag_view(a::AbstractVector) = a
@inline _diag_view(a::LinearAlgebra.Diagonal) = a.diag

# Build a callable `x -> a_shape(x)` for the given system under convention B, where
# `a_shape(x) = a(x) / s` and `s = L_1(a(x_0))/D` is the noise-magnitude scale extracted
# at the reference state `x_0 = current_state(sys)`. Returns:
#
#   * `(nothing, :additive)`            for additive noise (a is constant; callers use
#                                       their own precomputed `A = inv(a_shape)`).
#   * `(callable -> Vector, :diagonal)` when `σσ'` is diagonal at the reference state.
#   * `(callable -> Matrix, :general)`  otherwise.
function _normalized_a_shape_func(sys::CoupledSDEs)
    sys.noise_type[:additive] && return nothing, :additive
    g = diffusion_function(sys)
    params = current_parameters(sys)
    σ0 = g(current_state(sys), params, 0.0)
    a0 = σ0 * σ0'
    s = LinearAlgebra.norm(a0, 1) / size(a0, 1)
    if LinearAlgebra.isdiag(a0)
        f = let g = g, params = params, s = s
            x -> begin
                σ = g(x, params, 0.0)
                M = σ * σ'
                [M[k, k] / s for k in 1:size(M, 1)]
            end
        end
        return f, :diagonal
    else
        f = let g = g, params = params, s = s
            x -> begin
                σ = g(x, params, 0.0)
                Matrix((σ * σ') / s)
            end
        end
        return f, :general
    end
end

function path_velocity!(v, path, time; order = 4)
    if order == 2
        @views begin
            v[:, 1] .= (path[:, 2] .- path[:, 1]) / (time[2] - time[1])
            v[:, end] .= (path[:, end] .- path[:, end - 1]) / (time[end] - time[end - 1])
            for i in 2:(size(path, 2) - 1)
                v[:, i] .= (path[:, i + 1] .- path[:, i - 1]) / (time[i + 1] - time[i - 1])
            end
        end
    elseif order == 4
        @views begin
            v[:, 1] .= (path[:, 2] .- path[:, 1]) / (time[2] - time[1])
            v[:, end] .= (path[:, end] .- path[:, end - 1]) / (time[end] - time[end - 1])
            v[:, 2] .= (path[:, 3] .- path[:, 1]) / (time[3] - time[1])
            v[:, end - 1] .=
                (path[:, end] .- path[:, end - 2]) / (time[end] - time[end - 2])
            for i in 3:(size(path, 2) - 2)
                v[:, i] .= (
                    (
                        -path[:, i + 2] .+ 8 * path[:, i + 1] .- 8 * path[:, i - 1] .+
                            path[:, i - 2]
                    ) / (6 * (time[i + 1] - time[i - 1]))
                )
            end
        end
    end
    return v
end
