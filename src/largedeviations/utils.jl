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

function proper_MAM_system(ds::CoupledSDEs)
    for trait in (:additive, :invertible, :autonomous)
        if !ds.noise_type[trait]
            throw(
                ArgumentError(
                    "The minimal action method is only applicable for autonomous invertible additive noise. The noise type of the system is not $trait.",
                ),
            )
        end
    end
    return
end

"""
    proper_FW_system(ds::CoupledSDEs)

Validates that `ds` is a valid input for the Freidlin-Wentzell Hamiltonian path
(gMAM or sgMAM via `FreidlinWentzellHamiltonian(ds)`). Only requires autonomous noise;
rank-deficiency is detected by `_validate_and_classify_a` (called at workspace
construction with the path's reference state). Returns `nothing` on success.
"""
function proper_FW_system(ds::CoupledSDEs)
    if !ds.noise_type[:autonomous]
        throw(
            ArgumentError(
                "Freidlin-Wentzell methods are only applicable for autonomous noise.",
            ),
        )
    end
    return nothing
end

_fd_step(::Type{T}) where {T} = max(sqrt(eps(real(T))), real(T)(1.0e-8))

_as_diffusion_matrix(σ::AbstractMatrix) = σ
_as_diffusion_matrix(σ) = LinearAlgebra.Diagonal(σ)

const _RANK_DEFICIENT_MSG =
    "rank-deficient noise is not supported. Workarounds: add a small ε on the noiseless variable to make the covariance invertible, or supply a Hamiltonian directly via FreidlinWentzellHamiltonian{IIP, D}(H_x, H_p)."

function _check_rank!(a)
    M = a isa LinearAlgebra.Diagonal ? a : Matrix(a)
    if LinearAlgebra.cond(M) > 1 / sqrt(eps(real(eltype(M))))
        throw(ArgumentError(_RANK_DEFICIENT_MSG))
    end
    return nothing
end

_isdiag_numerical(::LinearAlgebra.Diagonal) = true
function _isdiag_numerical(a::AbstractMatrix; atol = sqrt(eps(real(eltype(a)))))
    @inbounds for j in axes(a, 2), i in axes(a, 1)
        i == j && continue
        abs(a[i, j]) > atol && return false
    end
    return true
end

"""
    _validate_and_classify_a(a, x_ref) -> is_diagonal::Bool

Validates `a(x)` at `x_ref` and at `x_ref ± h·e_l` for each coordinate `l`. Throws
`ArgumentError` if `a` is rank-deficient at any probe. Returns `true` if `a` is
numerically diagonal at every probe; `false` otherwise.
"""
function _validate_and_classify_a(a, x_ref::AbstractVector{T}) where {T}
    h = max(sqrt(eps(real(T))), real(T)(1.0e-6))
    D = length(x_ref)
    e = zeros(T, D)
    a0 = a(x_ref)
    _check_rank!(a0)
    is_diag = _isdiag_numerical(a0)
    @inbounds for l in 1:D
        e[l] = h
        ap = a(x_ref .+ e)
        _check_rank!(ap)
        is_diag &= _isdiag_numerical(ap)
        am = a(x_ref .- e)
        _check_rank!(am)
        is_diag &= _isdiag_numerical(am)
        e[l] = 0
    end
    return is_diag
end

"""
    _trace_normalized_a(ds::ContinuousTimeDynamicalSystem)

Returns a callable `x -> a(x)` where `a` is the trace-normalized diffusion tensor
`σ(x)σ(x)' / s` with `s = tr(σ(u₀)σ(u₀)')/D`. For constant noise, returns
`Base.Returns(a_const)` (constant). For `CoupledODEs`, returns `Returns(I)`.
"""
function _trace_normalized_a(ds::ContinuousTimeDynamicalSystem)
    D = dimension(ds)
    if ds isa CoupledODEs
        return Returns(LinearAlgebra.Diagonal(ones(Float64, D)))
    end
    σ_fn = diffusion_function(ds)
    ps = current_parameters(ds)
    σ0 = _as_diffusion_matrix(σ_fn(current_state(ds), ps, 0.0))
    a0 = σ0 * σ0'
    s = LinearAlgebra.tr(a0) / D
    if ds.noise_type[:additive]
        a_const = LinearAlgebra.isdiag(a0) ?
            LinearAlgebra.Diagonal(collect(LinearAlgebra.diag(a0)) ./ s) :
            Matrix(a0 ./ s)
        return Returns(a_const)
    end
    return let σ_fn = σ_fn, ps = ps, s = s
        x -> begin
            σx = _as_diffusion_matrix(σ_fn(x, ps, 0.0))
            (σx * σx') / s
        end
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

_thread_count() =
    hasmethod(Threads.nthreads, Tuple{Symbol}) ?
    (Threads.nthreads(:default) + Threads.nthreads(:interactive)) :
    Threads.nthreads()
