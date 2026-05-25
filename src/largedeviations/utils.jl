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
function interpolate_path!(path::Matrix, α, s, scratch::Vector = similar(α))
    D, N = size(path)
    @inbounds α[1] = zero(eltype(α))
    @inbounds for i in 2:N
        d2 = zero(eltype(α))
        for k in 1:D
            δ = path[k, i] - path[k, i - 1]
            d2 += δ * δ
        end
        α[i] = α[i - 1] + sqrt(d2)
    end
    if α[N] <= eps(real(eltype(α)))
        throw(
            ArgumentError(
                "interpolate_path!: path has zero arclength (all points coincide); cannot normalize.",
            ),
        )
    end
    invL = inv(α[N])
    @inbounds for i in 1:N
        α[i] *= invL
    end
    @inbounds α[N] = one(eltype(α))  # force exact endpoint to avoid round-off out-of-domain
    @inbounds for dof in 1:D
        linear_interp!(scratch, α, view(path, dof, :), s)
        for k in 1:N
            path[dof, k] = scratch[k]
        end
    end
    return nothing
end

"""
    proper_FW_system(ds::CoupledSDEs)

Validates that `ds` is a valid input for the Freidlin-Wentzell path methods (MAM via
[`minimize_action`](@ref), gMAM via [`minimize_geometric_action`](@ref), sgMAM via
[`FreidlinWentzellHamiltonian`](@ref) + `minimize_geometric_action`). Only requires
`:autonomous` noise; multiplicative and state-dependent diffusion are admissible.
Rank-deficiency of `a(x) = σ(x)σ(x)ᵀ` is detected separately by `_validate_and_classify_a`
(called at workspace / cache construction with the path's reference state). Returns
`nothing` on success; throws `ArgumentError` otherwise.
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

Validates `a(x)` at `x_ref` and at `x_ref ± h·eₗ` for each coordinate `l`, where `h` is a
finite-difference step. Throws `ArgumentError` if `a` is rank-deficient at any probe.
Returns `true` if `a` is numerically diagonal at every probe; `false` otherwise. The
return value drives `Val{true}/Val{false}` dispatch on the decoupled-vs-coupled cache
path in `build_sgmam_cache` and `geometric_gradient_workspace`.
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

Returns a callable `x -> a(x)` where `a(x) = σ(x) σ(x)ᵀ / s` is the trace-normalized
diffusion tensor and `s = tr(σ(u₀) σ(u₀)ᵀ) / D` is computed once at the current state `u₀
= current_state(ds)`. The normalization makes the Freidlin-Wentzell action invariant to
overall noise rescaling, so values returned by [`fw_action`](@ref), [`om_action`](@ref),
and [`geometric_action`](@ref) do not depend on `noise_strength` under the FW limit.

Return shape:
* `CoupledODEs` → `Returns(Diagonal(ones(D)))` (identity metric).
* Additive noise → `Returns(a_const)` with `a_const::Diagonal` if `σσᵀ` is diagonal,
  else `a_const::Matrix`. The closure is a `Base.Returns`, so callers can dispatch on it.
* State-dependent noise → a closure `x -> (σ(x) σ(x)ᵀ) / s` returning `Matrix` (or
  `Diagonal` when `σ` is supplied as a vector).
"""
function _trace_normalized_a(ds::ContinuousTimeDynamicalSystem)
    D = dimension(ds)
    if ds isa CoupledODEs
        T = eltype(current_state(ds))
        return Returns(LinearAlgebra.Diagonal(ones(T, D)))
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

"""
    path_velocity!(v, path, time; order = 4) -> v
    path_velocity(path, time; order = 4) -> v

Compute the time derivative ``\\dot\\phi(t)`` of a discrete path by finite differences.
`path_velocity!` writes into the preallocated buffer `v`; `path_velocity` allocates a new
matrix the same size and shape as `path`.

# Arguments
* `path`: `D × N` matrix with the path points in columns (`D` is the state dimension,
  `N` the number of time samples).
* `time`: length-`N` vector of monotonically increasing time points. Spacing need not be
  uniform; the stencils use the actual time differences. For the `order = 4` stencil to be
  accurate the spacing should be (approximately) uniform; nonuniform `time` is still
  handled but only with the formal accuracy of the corresponding uniform stencil.
* `v` (in-place form only): `D × N` matrix; will be overwritten in place.

# Keyword arguments
* `order::Int = 4`: order of the central finite-difference stencil used for interior
  points. Must be `2` or `4`. Other values are no-ops (the buffer is returned untouched).
  - `order = 2`: 3-point central differences at interior points; 1st-order forward /
    backward differences at the two endpoints.
  - `order = 4`: 5-point central differences for `i ∈ 3:N-2`; 2nd-order central
    differences at `i = 2` and `i = N-1`; 1st-order forward / backward differences at the
    endpoints `i = 1` and `i = N`.

# Returns
The velocity matrix `v` (same shape as `path`).
"""
function path_velocity!(v, path, time; order = 4)
    N = size(path, 2)
    if order == 2
        @inbounds @views begin
            inv_h1 = 1 / (time[2] - time[1])
            inv_hN = 1 / (time[end] - time[end - 1])
            @. v[:, 1] = (path[:, 2] - path[:, 1]) * inv_h1
            @. v[:, end] = (path[:, end] - path[:, end - 1]) * inv_hN
            for i in 2:(N - 1)
                inv_hi = 1 / (time[i + 1] - time[i - 1])
                @. v[:, i] = (path[:, i + 1] - path[:, i - 1]) * inv_hi
            end
        end
    elseif order == 4
        @inbounds @views begin
            inv_h1 = 1 / (time[2] - time[1])
            inv_hN = 1 / (time[end] - time[end - 1])
            inv_h2 = 1 / (time[3] - time[1])
            inv_hM = 1 / (time[end] - time[end - 2])
            @. v[:, 1] = (path[:, 2] - path[:, 1]) * inv_h1
            @. v[:, end] = (path[:, end] - path[:, end - 1]) * inv_hN
            @. v[:, 2] = (path[:, 3] - path[:, 1]) * inv_h2
            @. v[:, end - 1] = (path[:, end] - path[:, end - 2]) * inv_hM
            for i in 3:(N - 2)
                inv6 = 1 / (6 * (time[i + 1] - time[i - 1]))
                @. v[:, i] = (-path[:, i + 2] + 8 * path[:, i + 1] - 8 * path[:, i - 1] + path[:, i - 2]) * inv6
            end
        end
    end
    return v
end

"""
    path_velocity(path, time; order = 4) -> v

Allocating variant of [`path_velocity!`](@ref). Returns a freshly allocated `D × N` matrix
of velocities. See [`path_velocity!`](@ref) for argument semantics and the supported
stencils.
"""
path_velocity(path, time; order = 4) =
    path_velocity!(zeros(eltype(path), size(path)), path, time; order)
