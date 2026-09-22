"""
$(TYPEDEF)

Configuration for the Chebyshev Ritz discretization of a Noether-reduced on-shell action.

The path is represented by its values at `degree + 1` Chebyshev-Lobatto points. The
on-shell action is evaluated on a separate, oversampled Clenshaw-Curtis grid and then
minimized with an Optimization.jl-compatible optimizer.

For the Freidlin-Wentzell functional, `energy = 0` gives exactly the geometric action.
Nonzero `energy` selects a finite-time shell of the time-translation Noether invariant;
the shell-reality condition is checked while evaluating the action. The same reduction is
also available for the additive-noise Onsager-Machlup functional through the
`functional = "OM"` keyword of [`minimize_geometric_action`](@ref).

# Keyword arguments
  - `optimizer = Optimisers.Adam()`: Optimization.jl-compatible nonlinear optimizer.
  - `degree::Int = 8`: degree of the global Chebyshev path polynomial.
  - `quadrature_points::Int = 10 * degree`: Clenshaw-Curtis quadrature degree. The
    quadrature uses `quadrature_points + 1` nodes.
  - `energy::Real = 0.0`: finite on-shell Noether energy. Signed values are allowed when
    the corresponding shell remains real along the path.
  - `ad_type = OptimizationBase.AutoFiniteDiff()`: automatic-differentiation backend
    used by Optimization.jl.
"""
struct Ritz{O, A, T <: Real}
    optimizer::O
    degree::Int
    quadrature_points::Int
    energy::T
    ad_type::A
end

function Ritz(;
        optimizer = Optimisers.Adam(),
        degree::Int = 8,
        quadrature_points::Int = 10 * degree,
        energy::Real = 0.0,
        ad_type = OptimizationBase.AutoFiniteDiff(),
    )
    degree >= 2 || throw(ArgumentError("Ritz degree must be at least 2"))
    quadrature_points > degree ||
        throw(ArgumentError("quadrature_points must be greater than degree"))
    isfinite(energy) || throw(ArgumentError("Ritz energy must be finite"))
    E = float(energy)
    return Ritz{typeof(optimizer), typeof(ad_type), typeof(E)}(
        optimizer, degree, quadrature_points, E, ad_type
    )
end

# Increasing Chebyshev-Lobatto nodes on [-1, 1], as used in the Ritz paper.
function _ritz_chebyshev_nodes(n::Int, ::Type{T} = Float64) where {T <: AbstractFloat}
    return [-cospi(T(j) / T(n)) for j in 0:n]
end

# Barycentric weights for the increasing Chebyshev-Lobatto nodes above. A common
# multiplicative factor is irrelevant.
function _ritz_barycentric_weights(n::Int, ::Type{T} = Float64) where {T <: AbstractFloat}
    w = Vector{T}(undef, n + 1)
    @inbounds for j in 0:n
        endpoint = (j == 0 || j == n) ? T(0.5) : one(T)
        w[j + 1] = isodd(j) ? -endpoint : endpoint
    end
    return w
end

function _ritz_differentiation_matrix(nodes, bary_weights)
    N = length(nodes)
    T = promote_type(eltype(nodes), eltype(bary_weights))
    D = zeros(T, N, N)
    @inbounds for i in 1:N
        rowsum = zero(T)
        for j in 1:N
            i == j && continue
            Dij = bary_weights[j] / (bary_weights[i] * (nodes[i] - nodes[j]))
            D[i, j] = Dij
            rowsum += Dij
        end
        D[i, i] = -rowsum
    end
    return D
end

function _ritz_barycentric_matrix(nodes, bary_weights, eval_nodes)
    T = promote_type(eltype(nodes), eltype(bary_weights), eltype(eval_nodes))
    B = zeros(T, length(eval_nodes), length(nodes))
    tmp = Vector{T}(undef, length(nodes))
    tol = T(32) * eps(T)

    @inbounds for i in eachindex(eval_nodes)
        v = eval_nodes[i]
        matched = 0
        for j in eachindex(nodes)
            if abs(v - nodes[j]) <= tol * max(one(T), abs(v), abs(nodes[j]))
                matched = j
                break
            end
        end
        if matched != 0
            B[i, matched] = one(T)
            continue
        end

        denom = zero(T)
        for j in eachindex(nodes)
            tj = bary_weights[j] / (v - nodes[j])
            tmp[j] = tj
            denom += tj
        end
        for j in eachindex(nodes)
            B[i, j] = tmp[j] / denom
        end
    end
    return B
end

# Clenshaw-Curtis quadrature on [-1, 1], following Trefethen's Spectral Methods
# in MATLAB. The ordering is increasing to match _ritz_chebyshev_nodes.
function _ritz_clenshaw_curtis(n::Int, ::Type{T} = Float64) where {T <: AbstractFloat}
    n >= 2 || throw(ArgumentError("Clenshaw-Curtis degree must be at least 2"))
    nodes = _ritz_chebyshev_nodes(n, T)
    weights = zeros(T, n + 1)
    nT = T(n)

    if iseven(n)
        weights[1] = weights[end] = inv(T(n^2 - 1))
        @inbounds for j in 1:(n - 1)
            θ = T(pi) * T(j) / nT
            v = one(T)
            for k in 1:(n ÷ 2 - 1)
                v -= T(2) * cos(T(2 * k) * θ) / T(4 * k^2 - 1)
            end
            v -= cos(T(n) * θ) / T(n^2 - 1)
            weights[j + 1] = T(2) * v / nT
        end
    else
        weights[1] = weights[end] = inv(T(n^2))
        @inbounds for j in 1:(n - 1)
            θ = T(pi) * T(j) / nT
            v = one(T)
            for k in 1:((n - 1) ÷ 2)
                v -= T(2) * cos(T(2 * k) * θ) / T(4 * k^2 - 1)
            end
            weights[j + 1] = T(2) * v / nT
        end
    end
    return nodes, weights
end

function _ritz_spectral_matrices(degree::Int, quadrature_points::Int, ::Type{T}) where {T <: AbstractFloat}
    nodes = _ritz_chebyshev_nodes(degree, T)
    bary = _ritz_barycentric_weights(degree, T)
    D = _ritz_differentiation_matrix(nodes, bary)
    qnodes, qweights = _ritz_clenshaw_curtis(quadrature_points, T)
    B = _ritz_barycentric_matrix(nodes, bary, qnodes)
    C = B * D
    return nodes, bary, B, C, qweights
end

function _ritz_functional(functional)
    F = Symbol(functional)
    F in (:FW, :OM) || throw(
        ArgumentError(
            "Unknown Ritz action functional `$(F)`. Supported values are \"FW\" and \"OM\".",
        ),
    )
    return Val(F)
end

_ritz_validate_functional(::CoupledSDEs, ::Val{:FW}, _) = nothing
function _ritz_validate_functional(sys::CoupledSDEs, ::Val{:OM}, noise_strength)
    sys.noise_type[:additive] || throw(
        ArgumentError(
            "Ritz with functional = \"OM\" is implemented only for additive noise, matching om_action. Use functional = \"FW\" for state-dependent / multiplicative noise.",
        ),
    )
    noise_strength isa Real || throw(
        ArgumentError("noise_strength must be specified as a real number for functional = \"OM\""),
    )
    noise_strength >= 0 || throw(ArgumentError("noise_strength must be nonnegative"))
    return nothing
end

_ritz_action_auxiliary(::CoupledSDEs, ::Val{:FW}) = nothing
_ritz_action_auxiliary(sys::CoupledSDEs, ::Val{:OM}) = (jacobian(sys), sys.p0)

_ritz_velocity_independent_term(::Val{:FW}, _, _, _) = 0
function _ritz_velocity_independent_term(::Val{:OM}, aux, x, noise_strength)
    jac, p = aux
    σ = oftype(x[1], noise_strength)
    return σ^2 * tr(jac(x, p, 0)) / 2
end

# Backwards-compatible private FW entry point used by the spectral certification tests.
function _ritz_on_shell_action(sys, path, B, C, weights, energy, A_at)
    return _ritz_on_shell_action(
        sys, path, B, C, weights, energy, A_at, Val(:FW), nothing
    )
end

function _ritz_on_shell_action(
        sys, path, B, C, weights, energy, A_at, functional, noise_strength,
    )
    qpath = path * transpose(B)
    qvelocity = path * transpose(C)
    S = zero(eltype(qpath))
    E = oftype(S, energy)
    aux = _ritz_action_auxiliary(sys, functional)

    @views @inbounds for j in axes(qpath, 2)
        x = qpath[:, j]
        xp = qvelocity[:, j]
        b = drift(sys, x)
        A = _eval_metric(A_at, x)
        b2 = dot(b, A, b)
        xp2 = dot(xp, A, xp)
        cross = dot(xp, A, b)
        U = _ritz_velocity_independent_term(functional, aux, x, noise_strength)

        # For L = 1/2 |ẋ-b|_A^2 + U(x), time-translation invariance gives
        # E = 1/2 (|ẋ|_A^2 - |b|_A^2) - U. Eliminating dt/du therefore gives
        # |ẋ|_A^2 = |b|_A^2 + 2(E+U). U=0 is FW; for OM,
        # U = σ² div(b)/2 and may make the original Lagrangian non-positive.
        shell_speed2 = b2 + 2 * (E + U)
        shell_speed2 >= zero(shell_speed2) || throw(
            DomainError(
                shell_speed2,
                "Ritz on-shell action is not real on the requested energy shell; change energy or the path",
            ),
        )

        # At E = 0 the numerator equals shell_speed2, so writing the prefactor as
        # sqrt(shell_speed2) keeps zero-speed endpoints regular. For E != 0 a zero
        # shell speed is a genuine singularity of the Noether elimination.
        prefactor = if iszero(E)
            sqrt(shell_speed2)
        else
            iszero(shell_speed2) && throw(
                DomainError(
                    shell_speed2,
                    "nonzero-energy Ritz shell has zero speed and is singular",
                ),
            )
            (E + b2 + 2U) / sqrt(shell_speed2)
        end
        S += weights[j] * (prefactor * sqrt(xp2) - cross)
    end
    return S
end

function _ritz_pack_path(z, x_i, x_f, D::Int, degree::Int)
    interior = reshape(z, D, degree - 1)
    return hcat(x_i, interior, x_f)
end

# Resample an ordinary initial path, assumed uniformly parameterized on [-1, 1], onto
# the Chebyshev-Lobatto degrees of freedom. This affects only the optimizer initial guess.
function _ritz_resample_initial(init::AbstractMatrix, nodes)
    D, M = size(init)
    M >= 2 || throw(ArgumentError("initial path must contain at least two points"))
    T = promote_type(eltype(init), eltype(nodes))
    out = Matrix{T}(undef, D, length(nodes))
    scale = T(M - 1) / T(2)

    @inbounds for j in eachindex(nodes)
        s = (nodes[j] + one(T)) * scale
        left0 = clamp(floor(Int, s), 0, M - 1)
        if left0 == M - 1
            for d in 1:D
                out[d, j] = init[d, M]
            end
        else
            τ = s - T(left0)
            left = left0 + 1
            right = left + 1
            for d in 1:D
                out[d, j] = (one(T) - τ) * init[d, left] + τ * init[d, right]
            end
        end
    end
    return out
end

function _ritz_straight_initial(x_i, x_f, nodes, ::Type{T}) where {T <: AbstractFloat}
    D = length(x_i)
    path = Matrix{T}(undef, D, length(nodes))
    @inbounds for j in eachindex(nodes)
        τ = (nodes[j] + one(T)) / T(2)
        for d in 1:D
            path[d, j] = (one(T) - τ) * x_i[d] + τ * x_f[d]
        end
    end
    return path
end

function _minimize_geometric_action_ritz(
        sys::CoupledSDEs,
        init::AbstractMatrix,
        method::Ritz;
        functional = "FW",
        noise_strength = nothing,
        maxiters::Int = 1000,
        abstol::Real = NaN,
        reltol::Real = NaN,
        output_points::Int = 100,
        show_progress::Bool = true,
    )
    proper_FW_system(sys)
    output_points >= 2 || throw(ArgumentError("output_points must be at least 2"))
    size(init, 1) == dimension(sys) ||
        throw(DimensionMismatch("initial path dimension does not match system dimension"))
    functional = _ritz_functional(functional)
    _ritz_validate_functional(sys, functional, noise_strength)

    T = promote_type(Float64, eltype(init), typeof(method.energy))
    nodes, bary, B, C, qweights =
        _ritz_spectral_matrices(method.degree, method.quadrature_points, T)
    collocation_path = _ritz_resample_initial(init, nodes)
    x_i = copy(collocation_path[:, 1])
    x_f = copy(collocation_path[:, end])
    D = size(collocation_path, 1)

    _validate_and_classify_a(_trace_normalized_a(sys), collect(x_i))
    A_at = _action_metric(sys)
    z0 = vec(copy(collocation_path[:, 2:(end - 1)]))
    objective = let sys = sys, x_i = x_i, x_f = x_f, D = D,
            degree = method.degree, B = B, C = C, qweights = qweights,
            energy = method.energy, A_at = A_at, functional = functional,
            noise_strength = noise_strength
        (z, _) -> begin
            path = _ritz_pack_path(z, x_i, x_f, D, degree)
            _ritz_on_shell_action(
                sys, path, B, C, qweights, energy, A_at, functional, noise_strength
            )
        end
    end

    optf = SciMLBase.OptimizationFunction(objective, method.ad_type)
    prob = SciMLBase.OptimizationProblem(optf, z0, ())
    progress = Progress(maxiters; enabled = show_progress)
    function callback(_, _)
        show_progress && next!(progress)
        return false
    end
    sol = solve(prob, method.optimizer; maxiters, callback, abstol, reltol)

    collocation_path = _ritz_pack_path(sol.u, x_i, x_f, D, method.degree)
    output_nodes = collect(range(-one(T), one(T); length = output_points))
    Bout = _ritz_barycentric_matrix(nodes, bary, output_nodes)
    output_path = collocation_path * transpose(Bout)
    S = objective(sol.u, nothing)
    return MinimumActionPath(StateSpaceSet(output_path'), S)
end

"""
    minimize_geometric_action(sys::CoupledSDEs, x_i, x_f, method::Ritz; kwargs...)
    minimize_geometric_action(sys::CoupledSDEs, init, method::Ritz; kwargs...)

Minimize a Noether-reduced on-shell action with a global Chebyshev Ritz method.

Unlike gMAM, the Ritz method does not evolve and repeatedly reparameterize a local path
grid. It represents the complete path by a degree-`n` Chebyshev interpolant, evaluates the
action by oversampled Clenshaw-Curtis quadrature, and minimizes directly over the interior
Chebyshev values. Endpoint constraints are therefore exact and the optimization itself is
unconstrained.

`functional = "FW"` (the default) selects the Freidlin-Wentzell functional. At
`energy = 0` its on-shell form is exactly the geometric Freidlin-Wentzell action; nonzero
`energy` selects a finite-time Noether shell. Signed energies are allowed whenever the
shell condition `|b|_Q² + 2E >= 0` remains satisfied along the path.

`functional = "OM"` selects the additive-noise Onsager-Machlup functional already exposed
by [`om_action`](@ref). Supply the physical `noise_strength` explicitly. Its
velocity-independent term `σ² div(b)/2` is retained in the Noether reduction, so the Ritz
formulation remains valid even where the Onsager-Machlup Lagrangian itself is not positive
definite. The requested shell must remain real,
`|b|_Q² + 2E + σ² div(b) >= 0`, along the path; otherwise a `DomainError` is thrown.

When an initial `D × N` path is supplied, it is treated as a uniformly parameterized
initial guess and resampled onto the Chebyshev-Lobatto points. The returned
[`MinimumActionPath`](@ref) is sampled uniformly in the path parameter using
`output_points` points (default 100).

See Kikuchi, Singh, Cates & Adhikari, *Phys. Rev. Research* **2**, 033208 (2020),
DOI 10.1103/PhysRevResearch.2.033208 for the Chebyshev Ritz discretization and the
Freidlin-Wentzell on-shell construction.

# Keyword arguments
  - `functional = "FW"`: action functional, either `"FW"` or `"OM"`.
  - `noise_strength = nothing`: required for `functional = "OM"`; ignored for `"FW"`.
  - `maxiters = 1000`: maximum nonlinear-optimization iterations.
  - `abstol = NaN`, `reltol = NaN`: optimization stopping tolerances.
  - `output_points = 100`: number of uniformly parameterized output path points.
  - `show_progress = true`: display optimization progress.
"""
function minimize_geometric_action(
        sys::CoupledSDEs, x_i, x_f, method::Ritz; kwargs...
    )
    length(x_i) == dimension(sys) ||
        throw(DimensionMismatch("x_i dimension does not match system dimension"))
    length(x_f) == dimension(sys) ||
        throw(DimensionMismatch("x_f dimension does not match system dimension"))
    T = promote_type(Float64, eltype(x_i), eltype(x_f), typeof(method.energy))
    nodes = _ritz_chebyshev_nodes(method.degree, T)
    init = _ritz_straight_initial(x_i, x_f, nodes, T)
    return _minimize_geometric_action_ritz(sys, init, method; kwargs...)
end

function minimize_geometric_action(
        sys::CoupledSDEs, init::AbstractMatrix, method::Ritz; kwargs...
    )
    return _minimize_geometric_action_ritz(sys, Matrix(init), method; kwargs...)
end

function minimize_geometric_action(
        sys::CoupledSDEs, init::StateSpaceSets.AbstractStateSpaceSet, method::Ritz; kwargs...
    )
    return minimize_geometric_action(sys, _path_matrix(init), method; kwargs...)
end
