"""
    _action_metric(sys::CoupledSDEs)

Returns the inverse of the trace-normalized diffusion tensor.
For constant `a`, returns the inverse matrix directly. For state-dependent `a`,
returns a callable `x -> inv(a(x))`. Callers route through `_eval_metric`.
"""
_action_metric(sys::CoupledSDEs) = _action_metric(_trace_normalized_a(sys), sys)
_action_metric(a::Base.Returns, sys::CoupledSDEs) = inv(a(current_state(sys)))
_action_metric(a, sys::CoupledSDEs) = let a = a
    x -> inv(a(x))
end

# `A` may be a matrix (constant metric) or callable `x -> A(x)` (state-dep).
_eval_metric(A::AbstractMatrix, _) = A
_eval_metric(A, x) = A(x)

"""
$(TYPEDSIGNATURES)

Freidlin-Wentzell action of `path` (a `D × N` matrix with `D = dimension(sys)`) with time
points `time` (length `N`) for the drift `b = dynamic_rule(sys)`:

```math
S_T[\\phi] \\;=\\; \\tfrac{1}{2}\\int_0^T \\big\\| \\dot\\phi - \\mathbf{b}(\\phi)\\big\\|_{\\mathbf{Q}(\\phi)}^2 \\,\\mathrm{d}t
```

where ``\\|a\\|_\\mathbf{Q}^2 = \\langle a, \\mathbf{Q}^{-1} a\\rangle`` (see `anorm`), ``T`` is the
total time of the path, and ``\\mathbf{Q}(x)`` is the trace-normalized diffusion tensor
`a(x) / s` with `s = tr(a(u₀))/D` pinned at `current_state(ds)`; for state-dependent
diffusions `s` depends on the chosen `u₀`.
"""
function fw_action(sys::CoupledSDEs, path, time)
    @assert all(diff(time) .≈ diff(time[1:2])) "Freidlin-Wentzell action is only defined for equispaced time"
    A_at = _action_metric(sys)
    integrand = fw_integrand(sys, path, time, A_at)
    S = zero(eltype(path))
    @inbounds for i in 1:(size(path, 2) - 1)
        S += (integrand[i + 1] + integrand[i]) / 2 * (time[i + 1] - time[i])
    end
    return S / 2
end

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
    if !sys.noise_type[:additive]
        throw(
            ArgumentError(
                "om_action currently implements the constant-diffusion Onsager-Machlup correction term and is only defined for additive noise. Use fw_action for the leading-order rate function under state-dependent / multiplicative noise.",
            ),
        )
    end
    σ = noise_strength
    S = zero(eltype(path))
    @views @inbounds for i in 1:(size(path, 2) - 1)
        S += σ^2 / 2 * (
            (div_drift(sys, path[:, i + 1]) + div_drift(sys, path[:, i])) / 2 *
                (time[i + 1] - time[i])
        )
    end
    return fw_action(sys, path, time) + S
end

"""
$(TYPEDSIGNATURES)

Computes the action functional specified by `functional` for a given CoupledSDEs `sys` and
`path` parameterized by `time`.

* `functional = "FW"`: Returns the Freidlin-Wentzell action ([`fw_action`](@ref))
* `functional = "OM"`: Returns the Onsager-Machlup action ([`om_action`](@ref))
"""
function action(sys::CoupledSDEs, path::Matrix, time, functional::AbstractString; noise_strength = nothing)
    return action(sys, path, time, Val{Symbol(functional)}(); noise_strength)
end

action(sys::CoupledSDEs, path::Matrix, time, ::Val{:FW}; noise_strength = nothing) =
    fw_action(sys, path, time)
action(sys::CoupledSDEs, path::Matrix, time, ::Val{:OM}; noise_strength = nothing) =
    om_action(sys, path, time, noise_strength)
function action(sys::CoupledSDEs, path::Matrix, time, ::Val{S}; noise_strength = nothing) where {S}
    throw(
        ArgumentError(
            "Unknown action functional `$(S)`. Supported values are \"FW\" and \"OM\".",
        ),
    )
end

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

function _geometric_action_from_drift(b::Function, path, arclength::Real, A)
    N = size(path, 2)
    T = eltype(path)
    v_buf = similar(path)
    integrand_buf = zeros(T, N)
    return _geometric_action_from_drift!(b, path, arclength, A, v_buf, integrand_buf)
end

function _geometric_action_from_drift!(b::Function, path, arclength::Real, A, v_buf, integrand_buf)
    N = size(path, 2)
    T = eltype(path)
    path_velocity!(v_buf, path, range(zero(T), T(arclength); length = N); order = 4)
    @views @inbounds for i in 1:N
        xi = path[:, i]
        vi = v_buf[:, i]
        drift_i = b(xi)
        A_i = _eval_metric(A, xi)
        integrand_buf[i] = anorm(vi, A_i) * anorm(drift_i, A_i) - dot(vi, A_i, drift_i)
    end
    S = zero(eltype(path))
    @inbounds for i in 1:(N - 1)
        S += (integrand_buf[i + 1] + integrand_buf[i]) / 2
    end
    return S * arclength / (N - 1)
end

"""
$(TYPEDSIGNATURES)

Computes the squared ``A``-norm ``\\|\\dot\\phi_t - b\\|^2_A`` (see [`fw_action`](@ref)).
Returns a vector of length `N` with one value per time point.
"""
function fw_integrand(sys::CoupledSDEs, path, time, A)
    v = path_velocity(path, time; order = 4)
    sqnorm = zeros(eltype(path), size(path, 2))
    b(x) = drift(sys, x)
    diff_buf = similar(path, size(path, 1))
    @views @inbounds for i in axes(path, 2)
        xi = path[:, i]
        drift_i = b(xi)
        A_i = _eval_metric(A, xi)
        for k in eachindex(diff_buf)
            diff_buf[k] = v[k, i] - drift_i[k]
        end
        sqnorm[i] = anorm(diff_buf, A_i; square = true)
    end
    return sqnorm
end
