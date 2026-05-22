"""
    _action_metric(sys::CoupledSDEs)

Returns `x -> A(x)`, where `A(x)` is the inverse of the trace-normalized diffusion
tensor `a(x)/s` with `s = tr(a(u₀))/D`. For additive `sys` the callable is `Returns(A)`
(constant); for state-dependent it evaluates `σ(x)σ(x)ᵀ / s` per call. Throws
`ArgumentError` (via `_classify_noise_shape`) on rank-deficient or non-autonomous noise.
"""
function _action_metric(sys::CoupledSDEs)
    _classify_noise_shape(sys)
    σ_fn = diffusion_function(sys)
    u₀ = current_state(sys); ps = current_parameters(sys)
    a_of(x) = let σx = σ_fn(x, ps, 0.0)
        σ_mat = σx isa AbstractMatrix ? σx : LinearAlgebra.Diagonal(σx)
        σ_mat * σ_mat'
    end
    a0 = a_of(u₀)
    s = LinearAlgebra.tr(a0) / size(a0, 1)
    if sys.noise_type[:additive]
        return Returns(inv(a0 / s))
    else
        return x -> inv(a_of(x) / s)
    end
end

"""
$(TYPEDSIGNATURES)

Freidlin-Wentzell action of `path` (a `D × N` matrix with `D = dimension(sys)`) with time
points `time` (length `N`) for the drift `b = dynamic_rule(sys)`:

```math
S_T[\\phi] \\;=\\; \\tfrac{1}{2}\\int_0^T \\big\\| \\dot\\phi - \\mathbf{b}(\\phi)\\big\\|_{\\mathbf{Q}(\\phi)}^2 \\,\\mathrm{d}t
```

where ``\\|a\\|_\\mathbf{Q}^2 = \\langle a, \\mathbf{Q}^{-1} a\\rangle`` (see `anorm`), ``T`` is the
total time of the path, and ``\\mathbf{Q}(x)`` is the trace-normalized diffusion tensor
`a(x) / (tr a(u₀) / D)`. The convention makes the returned value independent of the
`noise_strength` keyword and invariant under orthogonal changes of basis, and supports
both state-independent (additive) and state-dependent noise.
"""
function fw_action(sys::CoupledSDEs, path, time)
    @assert all(diff(time) .≈ diff(time[1:2])) "Freidlin-Wentzell action is only defined for equispaced time"
    A_at = _action_metric(sys)
    integrand = fw_integrand(sys, path, time, A_at)
    S = 0
    for i in 1:(size(path, 2) - 1)
        S += (integrand[i + 1] + integrand[i]) / 2 * (time[i + 1] - time[i])
    end
    return S / 2
end;

"""
    om_action(sys::CoupledSDEs, path, time, noise_strength)

Onsager-Machlup action of `path` (a `D × N` matrix) with time points `time` (length `N`)
for the drift `b = dynamic_rule(sys)`, at the given `noise_strength` σ:

```math
S^\\sigma_T[\\phi] \\;=\\; \\tfrac{1}{2}\\int_0^T \\Big(\\big\\|\\dot\\phi - \\mathbf{b}(\\phi)\\big\\|_\\mathbf{Q}^2
    \\;+\\; \\sigma^2 \\,\\nabla\\cdot \\mathbf{b}(\\phi)\\Big)\\,\\mathrm{d}t
```

where ``\\|a\\|_\\mathbf{Q}^2 = \\langle a, \\mathbf{Q}^{-1} a\\rangle`` (see `anorm`), ``T`` is the
total time, and ``\\mathbf{Q}`` is the trace-normalized `covariance_matrix(sys)`. The first
term is exactly [`fw_action`](@ref) and is independent of the `noise_strength` keyword
under the trace-normalize convention. The second term is the Onsager-Machlup finite-noise
correction parameterized by the explicit `noise_strength` argument: pass the σ at which
you want ``-\\log P[\\phi]`` evaluated (typically the same `noise_strength` you used at
construction). As ``\\sigma \\to 0``, `om_action` → `fw_action`.
"""
function om_action(sys::CoupledSDEs, path, time, noise_strength)
    @assert all(diff(time) .≈ diff(time[1:2])) "Onsager-Machlup action is only defined for equispaced time"

    σ = noise_strength
    # Trapezoidal quadrature of (σ²/2) ∫ ∇·f dt
    S = 0.0
    for i in 1:(size(path, 2) - 1)
        S +=
            σ^2 / 2 * (
            (div_drift(sys, path[:, i + 1]) + div_drift(sys, path[:, i])) / 2 *
                (time[i + 1] - time[i])
        )
    end
    return fw_action(sys, path, time) + S
end;

"""
$(TYPEDSIGNATURES)

Computes the action functional specified by `functional` for a given CoupledSDEs `sys` and
`path` parameterized by `time`.

* `functional = "FW"`: Returns the Freidlin-Wentzell action ([`fw_action`](@ref))
* `functional = "OM"`: Returns the Onsager-Machlup action ([`om_action`](@ref))
"""
function action(sys::CoupledSDEs, path::Matrix, time, functional; noise_strength = nothing)
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

Geometric Freidlin-Wentzell action of `path` (a `D × N` matrix) with the given `arclength`,
for the drift `b = dynamic_rule(sys)`. Equals ``\\inf_T S_T[\\varphi]`` of the
time-parameterized [`fw_action`](@ref) on the instanton:

```math
\\bar S[\\varphi] \\;=\\; \\int_0^L \\Big( \\|\\varphi'\\|_\\mathbf{Q}\\,\\|\\mathbf{b}(\\varphi)\\|_\\mathbf{Q}
    \\;-\\; \\langle \\varphi', \\mathbf{b}(\\varphi)\\rangle_\\mathbf{Q} \\Big)\\,\\mathrm{d}s
```

where ``s`` is the arclength coordinate, ``L`` the arclength, and ``\\|a\\|_\\mathbf{Q}^2 =
\\langle a, b\\rangle_\\mathbf{Q} = a^\\top \\mathbf{Q}^{-1} b`` (see `anorm`). ``\\mathbf{Q}`` is the
trace-normalized `covariance_matrix(sys)`. As with [`fw_action`](@ref), the returned value
is independent of the `noise_strength` keyword and rotation-invariant.
"""
function geometric_action(sys::CoupledSDEs, path, arclength = 1.0)
    A_at = _action_metric(sys)
    b(x) = drift(sys, x)
    return _geometric_action_from_drift(b, path, arclength, A_at)
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
function geometric_action(b::Function, path, arclength = 1.0; A = nothing)
    if A === nothing
        D = size(path, 1)
        A = LinearAlgebra.Diagonal(ones(eltype(path), D))
    end
    return _geometric_action_from_drift(b, path, arclength, A)
end

# `A` may be either a matrix (constant metric) or a callable `x -> A(x)` (state-dep).
_eval_metric(A::AbstractMatrix, _) = A
_eval_metric(A, x) = A(x)

function _geometric_action_from_drift(b::Function, path, arclength::Real, A)
    N = size(path, 2)
    v = path_velocity(path, range(0, arclength; length = N); order = 4)

    integrand = zeros(eltype(path), N)
    for i in 1:N
        drift_i = b(path[:, i])
        A_i = _eval_metric(A, path[:, i])
        integrand[i] = anorm(v[:, i], A_i) * anorm(drift_i, A_i) - dot(v[:, i], A_i, drift_i)
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
    v = path_velocity(path, time; order = 4)
    sqnorm = zeros(size(path, 2))
    b(x) = drift(sys, x)
    for i in axes(path, 2)
        # assumes the drift is time independent
        drift_i = b(path[:, i])
        A_i = _eval_metric(A, path[:, i])
        sqnorm[i] = anorm(v[:, i] - drift_i, A_i; square = true)
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
function path_velocity(path, time; order = 4)
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
