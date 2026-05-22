"""
A structure representing an extended phase space system where your dissipative vector
field is embedded in a doubled dimensional phase space. The struct stores the partial
derivatives `H_x`, `H_p` of the Freidlin/Wentzell Hamiltonian

``H(x, p) = ⟨b(x), p⟩ + (1/2) ⟨p, a(x)·p⟩,``

a callable `a` for the diffusion tensor, and a `NoiseShape` type parameter `NS` that
encodes how the inner algorithm loops dispatch.
"""
struct FreidlinWentzellHamiltonian{IIP, D, Hx, Hp, AF, NS <: NoiseShape}
    H_x::Hx
    H_p::Hp
    a::AF
end

function FreidlinWentzellHamiltonian(ds::ContinuousTimeDynamicalSystem)
    D = dimension(ds)
    if ds isa CoupledSDEs
        proper_FW_system(ds)
    end
    NS = _classify_noise_shape(ds)

    σ_fn = ds isa CoupledSDEs ? diffusion_function(ds) : nothing
    ps   = current_parameters(ds)

    # Decision for `a(x)`:
    #  * No noise (CoupledODEs):                 a ≡ I, wrap in Returns.
    #  * Additive (state-independent) SDE:       a is constant; wrap the (trace-
    #    normalized) value in Returns regardless of diagonal-ness so update_p!/
    #    update_x! do not repeatedly re-evaluate σ_fn for a constant tensor.
    #  * State-dependent SDE:                    closure that evaluates σ(x)σ(x)ᵀ
    #    per call, trace-normalized at u₀.
    is_constant_a = (ds isa CoupledODEs) || (ds isa CoupledSDEs && ds.noise_type[:additive])

    a_callable = if ds isa CoupledODEs
        Returns(LinearAlgebra.Diagonal(ones(Float64, D)))
    elseif is_constant_a
        u₀ = current_state(ds)
        σ0 = σ_fn(u₀, ps, 0.0)
        σ_mat = σ0 isa AbstractMatrix ? σ0 : LinearAlgebra.Diagonal(σ0)
        a0 = σ_mat * σ_mat'
        s = LinearAlgebra.tr(a0) / D
        a_const = LinearAlgebra.isdiag(a0) ?
            LinearAlgebra.Diagonal(collect(LinearAlgebra.diag(a0)) ./ s) :
            Matrix(a0 ./ s)
        Returns(a_const)
    else
        u₀ = current_state(ds)
        σ0 = σ_fn(u₀, ps, 0.0)
        σ_mat0 = σ0 isa AbstractMatrix ? σ0 : LinearAlgebra.Diagonal(σ0)
        a0 = σ_mat0 * σ_mat0'
        s = LinearAlgebra.tr(a0) / D
        let σ_fn = σ_fn, ps = ps, s = s
            x -> begin
                σx = σ_fn(x, ps, 0.0)
                σ_mat = σx isa AbstractMatrix ? σx : LinearAlgebra.Diagonal(σx)
                (σ_mat * σ_mat') / s
            end
        end
    end

    f   = dynamic_rule(ds)
    jac = jacobian(ds)

    # Constant `a` (additive SDE or CoupledODEs) gives ∂ₓa ≡ 0, so skip the FD
    # section entirely; this also handles the constant non-diagonal Σ case
    # correctly even though it classifies as GeneralNoise.
    skip_ax_fd = is_constant_a

    function H_x(x, p)
        Hx = similar(x)
        Nt = size(x, 2); Dn = size(x, 1)
        h_fd = max(sqrt(eps(eltype(x))), eltype(x)(1e-8))
        for idx in 1:Nt
            xi = x[:, idx]
            pi_v = p[:, idx]
            jax = jac(xi, ps, 0.0)
            @inbounds for idc in 1:Dn
                Hx[idc, idx] = dot(jax[:, idc], pi_v)
            end
            if !skip_ax_fd
                e = zeros(eltype(xi), Dn)
                sample = a_callable(xi)
                is_diag = sample isa LinearAlgebra.Diagonal || LinearAlgebra.isdiag(sample)
                @inbounds for l in 1:Dn
                    fill!(e, 0); e[l] = h_fd
                    ap = a_callable(xi .+ e)
                    am = a_callable(xi .- e)
                    if is_diag
                        dla = (LinearAlgebra.diag(ap) .- LinearAlgebra.diag(am)) ./ (2 * h_fd)
                        Hx[l, idx] += 0.5 * dot(pi_v .^ 2, dla)
                    else
                        Hx[l, idx] += 0.5 * dot(pi_v, ((ap .- am) ./ (2 * h_fd)) * pi_v)
                    end
                end
            end
        end
        return Hx
    end

    function H_p(x, p)
        Hp = similar(x)
        for idx in 1:size(x, 2)
            a_x = a_callable(x[:, idx])
            Hp[:, idx] = a_x * p[:, idx] .+ f(x[:, idx], ps, 0.0)
        end
        return Hp
    end

    return FreidlinWentzellHamiltonian{
        isinplace(ds), D, typeof(H_x), typeof(H_p), typeof(a_callable), typeof(NS),
    }(H_x, H_p, a_callable)
end

function FreidlinWentzellHamiltonian{IIP, D}(
        H_x::Function, H_p::Function;
        a = Returns(LinearAlgebra.Diagonal(ones(Float64, D))),
    ) where {IIP, D}
    NS = _classify_user_a(a, D)
    return FreidlinWentzellHamiltonian{IIP, D, typeof(H_x), typeof(H_p), typeof(a), typeof(NS)}(
        H_x, H_p, a,
    )
end

function prettyprint(mlp::FreidlinWentzellHamiltonian{IIP, D, Hx, Hp, AF, NS}) where {IIP, D, Hx, Hp, AF, NS}
    return "Freidlin-Wentzell Hamiltonian on $D-dimensional state space ($(NS)) with $(IIP ? "in-place" : "out-of-place") H_x and H_p"
end

Base.show(io::IO, mlp::FreidlinWentzellHamiltonian) = print(io, prettyprint(mlp))

"""
$(TYPEDSIGNATURES)

Performs the simplified geometric Minimal Action Method (sgMAM) on the given system `sys`.
Our implementation is only valid for additive noise.

This method computes the optimal path in the phase space of a Hamiltonian system that
minimizes the Freidlin–Wentzell action. The Hamiltonian functions `H_x` and `H_p` define
the system's dynamics in a doubled phase. The initial state `x_initial` is evolved
iteratively using constrained gradient descent over a specified number of iterations. The
method can display a progress meter and will stop early if the absolute tolerance
`abstol` or relative tolerance `reltol` is achieved.

The function returns a [`MinimumActionPath`](@ref) containing the final path, the action value,
the Lagrange multipliers (`.λ`), the momentum (`.generalized_momentum`), and the state derivatives (`.path_velocity`).
The implementation is based on the work of [Grafke et al. (2019)](https://homepages.warwick.ac.uk/staff/T.Grafke/simplified-geometric-minimum-action-method-for-the-computation-of-instantons.html).

The optional positional argument `optimizer` controls step-size adaptation. It defaults to
`GeometricGradient(; stepsize=1e3)`, which enables backtracking step-size control with an
initial step of `1e3` (see [`GeometricGradient`](@ref)). Pass
`GeometricGradient(; max_backtracks=0)` to use a fixed step size.

The step size is configured via `GeometricGradient(; stepsize=...)`. When backtracking is
enabled, prefer a **large** initial step size: rejected steps are cheap and the controller
reduces the step size automatically, so starting large gives fast early progress without
sacrificing accuracy.

## Keyword arguments

  - `maxiters::Int=1000`: maximum number of *outer* iterations (path updates). When
    backtracking is enabled, each outer iteration may perform up to
    `optimizer.max_backtracks + 1` trial steps.
  - `show_progress::Bool=false`: if true, display a progress bar
  - `verbose::Bool=false`: if true, print additional output
  - `abstol::Real=NaN`: absolute tolerance for early stopping based on action change
  - `reltol::Real=NaN`: relative tolerance for early stopping based on action change
"""
function minimize_geometric_action(
        sys::FreidlinWentzellHamiltonian,
        x_initial::Matrix{T},
        optimizer::GeometricGradient = GeometricGradient(; stepsize = 1.0e3);
        maxiters::Int = 1000,
        show_progress::Bool = false,
        verbose::Bool = false,
        abstol::Real = NaN,
        reltol::Real = NaN,
    ) where {T}
    Nx, Nt = size(x_initial)
    s = range(0; stop = 1, length = Nt)
    x, p, pdot, xdot, lambda, alpha = init_allocation(x_initial, Nt)
    xdotdot = zeros(size(xdot))

    cache = _build_sgmam_cache(sys, eltype(x), Nx, Nt)

    x_prev = similar(x)

    # Ensure a consistent starting path for action comparisons
    interpolate_path!(x, alpha, s)
    _sgmam_refresh!(xdot, p, lambda, x, sys)
    initial_action = FW_action(xdot, p)

    function try_step!(ϵ)
        update!(x, xdot, xdotdot, p, pdot, lambda, sys, ϵ; cache)
        interpolate_path!(x, alpha, s)
        _sgmam_refresh!(xdot, p, lambda, x, sys)
        return FW_action(xdot, p)
    end
    save!() = copyto!(x_prev, x)
    function restore!()
        copyto!(x, x_prev)
        return _sgmam_refresh!(xdot, p, lambda, x, sys)
    end

    current_action, _ = backtracking_optimize!(
        optimizer,
        try_step!,
        save!,
        restore!,
        initial_action;
        maxiters,
        abstol,
        reltol,
        verbose,
        show_progress,
    )
    return MinimumActionPath(
        StateSpaceSet(x'),
        current_action;
        λ = lambda,
        generalized_momentum = p,
        path_velocity = xdot,
    )
end
function minimize_geometric_action(
        sys::FreidlinWentzellHamiltonian,
        x_initial::StateSpaceSet,
        optimizer::GMAMOptimizer = GeometricGradient(; stepsize = 1.0e3);
        kwargs...,
    )
    return minimize_geometric_action(
        sys, Matrix(Matrix(x_initial)'), optimizer; kwargs...
    )
end

"""
$(TYPEDSIGNATURES)

Adaptive multi-phase variant of the sgMAM projected-gradient method. See
[`AdaptiveGeometricGradient`](@ref) for the algorithm.

`maxiters` here counts the *effective* (kept) inner iterations the path has experienced
across all probe windows. Each probe window does `2 * probe_length` actual gradient
updates but advances the path by `probe_length` accepted iterations, so wall time is
roughly twice that of a fixed-step run with the same `maxiters`.
"""
function minimize_geometric_action(
        sys::FreidlinWentzellHamiltonian,
        x_initial::Matrix{T},
        optimizer::AdaptiveGeometricGradient;
        maxiters::Int = 1000,
        show_progress::Bool = false,
        verbose::Bool = false,
        abstol::Real = NaN,
        reltol::Real = NaN,
    ) where {T}
    Nx = size(x_initial, 1)
    Nt = size(x_initial, 2)
    s = range(0; stop = 1, length = Nt)
    x, p, pdot, xdot, lambda, alpha = init_allocation(x_initial, Nt)
    xdotdot = zeros(size(xdot))

    cache = _build_sgmam_cache(sys, eltype(x), Nx, Nt)

    x_start = similar(x)        # path snapshot at start of probe window
    x_big_result = similar(x)   # store big-probe result while running small probe

    interpolate_path!(x, alpha, s)
    _sgmam_refresh!(xdot, p, lambda, x, sys)
    Tϵ = typeof(optimizer.stepsize)
    current_action = Tϵ(FW_action(xdot, p))

    # Run `n` projected-gradient updates at fixed `ϵ`; return final action, or Inf
    # if any iteration produced a non-finite result. Closure captures all buffers.
    function run_probe!(ϵ, n)
        S = oftype(ϵ, NaN)
        for _ in 1:n
            update!(x, xdot, xdotdot, p, pdot, lambda, sys, ϵ; cache)
            interpolate_path!(x, alpha, s)
            _sgmam_refresh!(xdot, p, lambda, x, sys)
            S = oftype(ϵ, FW_action(xdot, p))
            isfinite(S) || return oftype(ϵ, Inf)
        end
        return S
    end

    stepsize = optimizer.stepsize
    probe_len = optimizer.probe_length

    iters_used = 0
    progress = Progress(maxiters; dt = 0.5, enabled = show_progress)

    while iters_used < maxiters
        S_prev = current_action
        n = min(probe_len, maxiters - iters_used)

        copyto!(x_start, x)
        ϵ_big = clamp(stepsize, optimizer.stepsize_min, optimizer.stepsize_max)
        S_big = run_probe!(ϵ_big, n)
        copyto!(x_big_result, x)

        copyto!(x, x_start)
        _sgmam_refresh!(xdot, p, lambda, x, sys)
        ϵ_small = clamp(
            stepsize * optimizer.shrink, optimizer.stepsize_min, optimizer.stepsize_max
        )
        S_small = run_probe!(ϵ_small, n)

        # A probe is usable only if its final action is finite and at most S_prev;
        # otherwise it ran away numerically and must be discarded.
        big_ok = isfinite(S_big) && S_big <= S_prev
        small_ok = isfinite(S_small) && S_small <= S_prev
        accepted = true
        if small_ok && (!big_ok || S_small < S_big)
            # small probe already left its result in x with consistent buffers
            current_action = S_small
            stepsize = max(optimizer.stepsize_min, stepsize * optimizer.shrink)
        elseif big_ok
            copyto!(x, x_big_result)
            _sgmam_refresh!(xdot, p, lambda, x, sys)
            current_action = S_big
            stepsize = min(optimizer.stepsize_max, stepsize * optimizer.grow)
        else
            copyto!(x, x_start)
            _sgmam_refresh!(xdot, p, lambda, x, sys)
            stepsize = max(optimizer.stepsize_min, stepsize * optimizer.shrink^2)
            accepted = false
            verbose &&
                @info "Probe rejected at iters_used=$iters_used (S_big=$S_big, S_small=$S_small, S_prev=$S_prev); shrinking stepsize to $stepsize."
            if stepsize <= optimizer.stepsize_min
                verbose && @info "stepsize hit stepsize_min; stopping."
                break
            end
        end

        iters_used += n

        # Skip the tolerance check on a rejected probe: current_action == S_prev would
        # otherwise trigger a false-positive convergence when the controller is in
        # trouble, not at a minimum.
        abs_change = accepted ? abs(current_action - S_prev) : oftype(current_action, Inf)
        rel_change = if accepted
            current_action == 0 ? abs_change : abs_change / abs(current_action)
        else
            oftype(current_action, Inf)
        end
        if (isfinite(abstol) && abs_change < abstol) ||
                (isfinite(reltol) && rel_change < reltol)
            verbose && @info "Converged after $iters_used iterations: abs=$abs_change, rel=$rel_change"
            break
        end

        next!(
            progress;
            step = n,
            showvalues = [
                ("iters_used", iters_used),
                ("action", round(current_action; sigdigits = 6)),
                ("stepsize", round(stepsize; sigdigits = 3)),
                ("Stol", round(rel_change; sigdigits = 3)),
            ],
        )
    end

    return MinimumActionPath(
        StateSpaceSet(x'),
        current_action;
        λ = lambda,
        generalized_momentum = p,
        path_velocity = xdot,
    )
end

function init_allocation(x_initial, Nt)
    # preallocate
    x = deepcopy(x_initial)
    p = zeros(size(x))
    pdot = zeros(size(x))
    xdot = zeros(size(x))
    lambda = zeros(1, Nt) # Lagrange multiplier
    alpha = zeros(Nt)
    return x, p, pdot, xdot, lambda, alpha
end

function update!(x, xdot, xdotdot, p, pdot, lambda, sys, ϵ; cache = nothing)
    # xdot, p, and lambda are assumed to be pre-computed and consistent with x
    # (via _sgmam_refresh! or central_diff! + update_p!)
    central_diff!(pdot, p)
    Hx = sys.H_x(x, p)
    central_diff!(xdotdot, xdot)
    return update_x!(x, lambda, pdot, xdotdot, Hx, sys, ϵ; cache)
end

function update_x!(x, λ, p′, x′′, Hx, sys::FreidlinWentzellHamiltonian, ϵ; cache = nothing)
    _update_x!(x, λ, p′, x′′, Hx, sys, ϵ, cache)
    return nothing
end

_sgmam_nthreads() =
    hasmethod(Threads.nthreads, Tuple{Symbol}) ?
        (Threads.nthreads(:default) + Threads.nthreads(:interactive)) :
        Threads.nthreads()

# Allocate `nthreads` LinearSolve caches with placeholder Tridiagonal storage of
# the right size. Each `_update_x!` call reuses these by mutating the aliased A
# and b in place and calling `LinearSolve.reinit!`.
function _sgmam_init_tridiag_caches(::Type{T}, L::Int, nthreads::Int) where {T}
    return map(1:nthreads) do _
        d0  = ones(T, L)
        du0 = zeros(T, L - 1)
        dl0 = zeros(T, L - 1)
        T0  = LinearAlgebra.Tridiagonal(dl0, d0, du0)
        rhs0 = zeros(T, L)
        init(
            LinearProblem(T0, rhs0), LUFactorization();
            alias = SciMLBase.LinearAliasSpecifier(; alias_A = true, alias_b = true),
        )
    end
end

function _update_x!(
        x, λ, p′, x′′, Hx,
        sys::FreidlinWentzellHamiltonian{IIP, D, Hx_t, Hp_t, AF, AdditiveNoise}, ϵ,
        _cache = nothing,
    ) where {IIP, D, Hx_t, Hp_t, AF}
    a_diag = LinearAlgebra.diag(sys.a(view(x, :, 1)))
    a_inv  = inv.(a_diag)

    Nx, Nt = size(x); xa = x[:, 1]; xb = x[:, end]; idxc = 2:(Nt - 1)
    L = Nt - 2

    nthreads = _sgmam_nthreads()
    caches = _sgmam_init_tridiag_caches(eltype(x), L, nthreads)

    Threads.@threads for dof in 1:Nx
        c = a_inv[dof]
        cache = caches[Threads.threadid()]
        Tmat = cache.A
        @inbounds for i in 1:L
            Tmat.d[i] = 1 + 2 * ϵ * λ[i + 1]^2 * c
        end
        @inbounds for i in 1:(L - 1)
            Tmat.du[i] = -ϵ * λ[i + 1]^2 * c
            Tmat.dl[i] = -ϵ * λ[i + 2]^2 * c
        end
        rhs = cache.b
        @inbounds for k in 1:L
            i = k + 1
            rhs[k] = x[dof, i] + ϵ * (λ[i] * p′[dof, i] + Hx[dof, i] - c * λ[i]^2 * x′′[dof, i])
        end
        rhs[1]   += ϵ * c * λ[2]^2 * xa[dof]
        rhs[end] += ϵ * c * λ[end - 1]^2 * xb[dof]
        LinearSolve.reinit!(cache; A = Tmat, b = rhs)
        solve!(cache)
        x[dof, idxc] .= cache.u
    end
    return nothing
end

function _update_x!(
        x, λ, p′, x′′, Hx,
        sys::FreidlinWentzellHamiltonian{IIP, D, Hx_t, Hp_t, AF, DiagonalNoise}, ϵ,
        _cache = nothing,
    ) where {IIP, D, Hx_t, Hp_t, AF}
    Nx, Nt = size(x)
    a_at = Matrix{eltype(x)}(undef, Nx, Nt)
    for t in 1:Nt
        a_at[:, t] .= LinearAlgebra.diag(sys.a(view(x, :, t)))
    end

    xa = x[:, 1]; xb = x[:, end]; idxc = 2:(Nt - 1)
    L = Nt - 2

    nthreads = _sgmam_nthreads()
    caches = _sgmam_init_tridiag_caches(eltype(x), L, nthreads)

    Threads.@threads for dof in 1:Nx
        cache = caches[Threads.threadid()]
        Tmat = cache.A
        @inbounds for i in 1:L
            ii = i + 1
            α_i = ϵ * λ[ii]^2 / a_at[dof, ii]
            Tmat.d[i] = 1 + 2 * α_i
        end
        @inbounds for i in 1:(L - 1)
            Tmat.du[i] = -ϵ * λ[i + 1]^2 / a_at[dof, i + 1]
            Tmat.dl[i] = -ϵ * λ[i + 2]^2 / a_at[dof, i + 2]
        end
        rhs = cache.b
        @inbounds for k in 1:L
            i = k + 1
            rhs[k] = x[dof, i] + ϵ * (
                λ[i] * p′[dof, i] + Hx[dof, i] - (λ[i]^2 / a_at[dof, i]) * x′′[dof, i]
            )
        end
        α_left  = ϵ * λ[2]^2 / a_at[dof, 2]
        α_right = ϵ * λ[Nt - 1]^2 / a_at[dof, Nt - 1]
        rhs[1]   += α_left  * xa[dof]
        rhs[end] += α_right * xb[dof]
        LinearSolve.reinit!(cache; A = Tmat, b = rhs)
        solve!(cache)
        x[dof, idxc] .= cache.u
    end
    return nothing
end

"""
Cross-iteration cache for the GeneralNoise sgMAM `_update_x!` path: holds the sparse
block-tridiagonal matrix `M`, an index map from logical (i_in, k1, k2) entries into
`M.nzval`, an RHS buffer, and a reusable `LinearSolve.LinearCache`. The sparsity
pattern is fixed (depends only on `Nt`, `Nx`); only `nzval` and `rhs` change per call.
"""
struct _SgMAMGeneralCache{T, LC}
    M::SparseArrays.SparseMatrixCSC{T, Int}
    diag_nzval_idx::Array{Int, 3}    # (Nx, Nx, N_in)
    left_nzval_idx::Matrix{Int}      # (Nx, N_in - 1): for i_in in 2:N_in
    right_nzval_idx::Matrix{Int}     # (Nx, N_in - 1): for i_in in 1:N_in-1
    rhs::Vector{T}
    linear_cache::LC
end

"""
Dispatcher that builds a cross-iteration cache for the inner `_update_x!` linear
solve when the noise shape benefits from one. Currently only `GeneralNoise` uses a
cache; other shapes return `nothing` (per-call buffers in their `_update_x!`).
"""
_build_sgmam_cache(::FreidlinWentzellHamiltonian, ::Type, ::Int, ::Int) = nothing
function _build_sgmam_cache(
        ::FreidlinWentzellHamiltonian{IIP, D, Hx, Hp, AF, GeneralNoise},
        ::Type{T}, Nx::Int, Nt::Int,
    ) where {IIP, D, Hx, Hp, AF, T}
    return _build_sgmam_general_cache(T, Nx, Nt)
end

function _find_nzval_idx(M::SparseArrays.SparseMatrixCSC, r::Int, c::Int)
    @inbounds for k in M.colptr[c]:(M.colptr[c + 1] - 1)
        if M.rowval[k] == r
            return k
        end
    end
    error("entry ($r, $c) is not in the sparsity pattern")
end

function _build_sgmam_general_cache(::Type{T}, Nx::Int, Nt::Int) where {T}
    N_in = Nt - 2
    n    = N_in * Nx
    nnz_est = N_in * Nx * Nx + 2 * (N_in - 1) * Nx

    Iv = Vector{Int}(undef, nnz_est)
    Jv = Vector{Int}(undef, nnz_est)
    Vv = ones(T, nnz_est)   # placeholder so M is non-singular at build time
    pos = 0
    @inbounds for i_in in 1:N_in
        rb = (i_in - 1) * Nx
        for k1 in 1:Nx, k2 in 1:Nx
            pos += 1; Iv[pos] = rb + k1; Jv[pos] = rb + k2
        end
        if i_in > 1
            base_L = (i_in - 2) * Nx
            for k in 1:Nx
                pos += 1; Iv[pos] = rb + k; Jv[pos] = base_L + k
            end
        end
        if i_in < N_in
            base_R = i_in * Nx
            for k in 1:Nx
                pos += 1; Iv[pos] = rb + k; Jv[pos] = base_R + k
            end
        end
    end
    M = SparseArrays.sparse(Iv, Jv, Vv, n, n)

    diag_idx  = Array{Int}(undef, Nx, Nx, N_in)
    left_idx  = Matrix{Int}(undef, Nx, N_in - 1)
    right_idx = Matrix{Int}(undef, Nx, N_in - 1)
    @inbounds for i_in in 1:N_in
        rb = (i_in - 1) * Nx
        for k1 in 1:Nx, k2 in 1:Nx
            diag_idx[k1, k2, i_in] = _find_nzval_idx(M, rb + k1, rb + k2)
        end
        if i_in > 1
            base_L = (i_in - 2) * Nx
            for k in 1:Nx
                left_idx[k, i_in - 1] = _find_nzval_idx(M, rb + k, base_L + k)
            end
        end
        if i_in < N_in
            base_R = i_in * Nx
            for k in 1:Nx
                right_idx[k, i_in] = _find_nzval_idx(M, rb + k, base_R + k)
            end
        end
    end

    rhs = zeros(T, n)
    fill!(M.nzval, zero(T))
    lc = init(LinearProblem(M, rhs), LUFactorization())
    return _SgMAMGeneralCache{T, typeof(lc)}(M, diag_idx, left_idx, right_idx, rhs, lc)
end

# Fallback dispatch when a GeneralNoise system is invoked without an explicit cache
# (rare path: direct `update_x!` call from tests). Builds a one-shot cache and uses it.
function _update_x!(
        x, λ, p′, x′′, Hx,
        sys::FreidlinWentzellHamiltonian{IIP, D, Hx_t, Hp_t, AF, GeneralNoise}, ϵ,
        ::Nothing = nothing,
    ) where {IIP, D, Hx_t, Hp_t, AF}
    Nx, Nt = size(x)
    cache = _build_sgmam_general_cache(eltype(x), Nx, Nt)
    return _update_x!(x, λ, p′, x′′, Hx, sys, ϵ, cache)
end

function _update_x!(
        x, λ, p′, x′′, Hx,
        sys::FreidlinWentzellHamiltonian{IIP, D, Hx_t, Hp_t, AF, GeneralNoise}, ϵ,
        cache::_SgMAMGeneralCache,
    ) where {IIP, D, Hx_t, Hp_t, AF}
    Nx, Nt = size(x)
    N_in = Nt - 2

    M = cache.M
    rhs = cache.rhs
    xa = view(x, :, 1); xb = view(x, :, Nt)

    @inbounds for i_in in 1:N_in
        i = i_in + 1
        rb = (i_in - 1) * Nx
        A_i = sys.a(view(x, :, i))
        λi2 = λ[i]^2

        for k1 in 1:Nx, k2 in 1:Nx
            v = A_i[k1, k2] + (k1 == k2 ? 2 * ϵ * λi2 : zero(eltype(x)))
            M.nzval[cache.diag_nzval_idx[k1, k2, i_in]] = v
        end
        if i_in > 1
            for k in 1:Nx
                M.nzval[cache.left_nzval_idx[k, i_in - 1]] = -ϵ * λi2
            end
        end
        if i_in < N_in
            for k in 1:Nx
                M.nzval[cache.right_nzval_idx[k, i_in]] = -ϵ * λi2
            end
        end

        for k in 1:Nx
            acc = zero(eltype(x))
            for k2 in 1:Nx
                acc += A_i[k, k2] * (x[k2, i] + ϵ * (λ[i] * p′[k2, i] + Hx[k2, i]))
            end
            v = acc - ϵ * λi2 * x′′[k, i]
            if i_in == 1
                v += ϵ * λi2 * xa[k]
            end
            if i_in == N_in
                v += ϵ * λi2 * xb[k]
            end
            rhs[rb + k] = v
        end
    end

    lc = cache.linear_cache
    LinearSolve.reinit!(lc; A = M, b = rhs)
    solve!(lc)
    sol = lc.u
    @inbounds for i_in in 1:N_in
        rb = (i_in - 1) * Nx
        for k in 1:Nx
            x[k, i_in + 1] = sol[rb + k]
        end
    end
    return nothing
end

function update_p!(p, lambda, x, xdot, sys::FreidlinWentzellHamiltonian)
    _update_p!(p, lambda, x, xdot, sys)
    return nothing
end

function _update_p!(
        p, lambda, x, xdot,
        sys::FreidlinWentzellHamiltonian{IIP, D, Hx, Hp, AF, AdditiveNoise},
    ) where {IIP, D, Hx, Hp, AF}
    b_ = sys.H_p(x, zero(x))
    a_diag = LinearAlgebra.diag(sys.a(view(x, :, 1)))
    a_inv  = inv.(a_diag)
    num = sum(b_ .^ 2 .* a_inv; dims = 1)
    den = sum(xdot .^ 2 .* a_inv; dims = 1)
    lambda .= sqrt.(num ./ den)
    lambda[1] = 0
    lambda[end] = 0
    p .= (lambda .* xdot .- b_) .* a_inv
    return nothing
end

function _update_p!(
        p, lambda, x, xdot,
        sys::FreidlinWentzellHamiltonian{IIP, D, Hx, Hp, AF, DiagonalNoise},
    ) where {IIP, D, Hx, Hp, AF}
    b_ = sys.H_p(x, zero(x))
    Nt = size(x, 2)
    for t in 1:Nt
        a_t   = sys.a(view(x, :, t))
        diag_t = LinearAlgebra.diag(a_t)
        inv_t  = inv.(diag_t)
        num_t = sum(b_[:, t] .^ 2 .* inv_t)
        den_t = sum(xdot[:, t] .^ 2 .* inv_t)
        λ_t = den_t > 1e-28 ? sqrt(num_t / den_t) : zero(eltype(b_))
        lambda[1, t] = isfinite(λ_t) ? λ_t : zero(eltype(b_))
        @inbounds for i in 1:D
            p[i, t] = (lambda[1, t] * xdot[i, t] - b_[i, t]) * inv_t[i]
        end
    end
    lambda[1, 1] = 0
    lambda[1, end] = 0
    return nothing
end

function _update_p!(
        p, lambda, x, xdot,
        sys::FreidlinWentzellHamiltonian{IIP, D, Hx, Hp, AF, GeneralNoise},
    ) where {IIP, D, Hx, Hp, AF}
    b_ = sys.H_p(x, zero(x))
    Nt = size(x, 2)
    inv_b    = Vector{eltype(b_)}(undef, D)
    inv_xdot = Vector{eltype(b_)}(undef, D)
    rhs      = Vector{eltype(b_)}(undef, D)
    for t in 1:Nt
        a_t = sys.a(view(x, :, t))
        inv_b    .= a_t \ Vector(view(b_, :, t))
        inv_xdot .= a_t \ Vector(view(xdot, :, t))
        num_t = dot(view(b_, :, t),   inv_b)
        den_t = dot(view(xdot, :, t), inv_xdot)
        λ_t = den_t > 1e-28 ? sqrt(num_t / den_t) : zero(eltype(b_))
        lambda[1, t] = isfinite(λ_t) ? λ_t : zero(eltype(b_))
        @inbounds for i in 1:D
            rhs[i] = lambda[1, t] * xdot[i, t] - b_[i, t]
        end
        p[:, t] .= a_t \ rhs
    end
    lambda[1, 1] = 0
    lambda[1, end] = 0
    return nothing
end

function central_diff!(xdot, x)
    # ̇xₙ = 0.5(xₙ₊₁ - xₙ₋₁) central finite difference
    @views xdot[:, 2:(end - 1)] .= 0.5 .* (x[:, 3:end] .- x[:, 1:(end - 2)])
    return nothing
end

function _sgmam_refresh!(xdot, p, lambda, x, sys)
    central_diff!(xdot, x)
    update_p!(p, lambda, x, xdot, sys)
    return nothing
end

# Freidlin-Wentzell action of an instanton path discretized over a uniform arclength
# parameter `s ∈ [0, 1]`. With `xdot` produced by `central_diff!` (which does not
# divide by Δs), `xdot[:, i] ≈ Δs · dϕ/ds`, so `dot(xdot, p) ≈ ∫ p · dϕ`. On the
# zero-energy shell `H = p·f + |p|²/2 = 0`, the identity `p · dϕ = (|p|²/2) dt`
# holds along the instanton, so `∫ p · dϕ = ½ ∫ |p|² dt = S_FW`. No additional `/2`.
FW_action(xdot, p) = dot(xdot, p)
