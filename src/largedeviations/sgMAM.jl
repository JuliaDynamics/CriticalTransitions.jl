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
function minimize_simple_geometric_action(
        sys::FreidlinWentzellHamiltonian,
        x_initial::Matrix{T},
        optimizer::GeometricGradient = GeometricGradient(; stepsize = 1.0e3);
        maxiters::Int = 1000,
        show_progress::Bool = false,
        verbose::Bool = false,
        abstol::Real = NaN,
        reltol::Real = NaN,
    ) where {T}
    H_p, H_x = sys.H_p, sys.H_x

    Nx, Nt = size(x_initial)
    s = range(0; stop = 1, length = Nt)
    x, p, pdot, xdot, lambda, alpha = init_allocation(x_initial, Nt)
    xdotdot = zeros(size(xdot))

    x_prev = similar(x)

    # Ensure a consistent starting path for action comparisons
    interpolate_path!(x, alpha, s)
    _sgmam_refresh!(xdot, p, lambda, x, H_p)
    initial_action = FW_action(xdot, p)

    function try_step!(ϵ)
        update!(x, xdot, xdotdot, p, pdot, lambda, H_x, H_p, ϵ)
        interpolate_path!(x, alpha, s)
        _sgmam_refresh!(xdot, p, lambda, x, H_p)
        return FW_action(xdot, p)
    end
    save!() = copyto!(x_prev, x)
    function restore!()
        copyto!(x, x_prev)
        return _sgmam_refresh!(xdot, p, lambda, x, H_p)
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
function minimize_simple_geometric_action(
        sys,
        x_initial::StateSpaceSet,
        optimizer::GMAMOptimizer = GeometricGradient(; stepsize = 1.0e3);
        kwargs...,
    )
    return minimize_simple_geometric_action(
        sys, Matrix(Matrix(x_initial)'), optimizer; kwargs...
    )
end
function minimize_simple_geometric_action(
        sys::ContinuousTimeDynamicalSystem,
        x_initial::Matrix{<:Real},
        optimizer::GMAMOptimizer = GeometricGradient(; stepsize = 1.0e3);
        kwargs...,
    )
    return minimize_simple_geometric_action(
        FreidlinWentzellHamiltonian(sys), Matrix(Matrix(x_initial)'), optimizer; kwargs...
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
function minimize_simple_geometric_action(
        sys::FreidlinWentzellHamiltonian,
        x_initial::Matrix{T},
        optimizer::AdaptiveGeometricGradient;
        maxiters::Int = 1000,
        show_progress::Bool = false,
        verbose::Bool = false,
        abstol::Real = NaN,
        reltol::Real = NaN,
    ) where {T}
    H_p, H_x = sys.H_p, sys.H_x

    Nt = size(x_initial, 2)
    s = range(0; stop = 1, length = Nt)
    x, p, pdot, xdot, lambda, alpha = init_allocation(x_initial, Nt)
    xdotdot = zeros(size(xdot))

    x_start = similar(x)        # path snapshot at start of probe window
    x_big_result = similar(x)   # store big-probe result while running small probe

    interpolate_path!(x, alpha, s)
    _sgmam_refresh!(xdot, p, lambda, x, H_p)
    Tϵ = typeof(optimizer.stepsize)
    current_action = Tϵ(FW_action(xdot, p))

    # Run `n` projected-gradient updates at fixed `ϵ`; return final action, or Inf
    # if any iteration produced a non-finite result. Closure captures all buffers.
    function run_probe!(ϵ, n)
        S = oftype(ϵ, NaN)
        for _ in 1:n
            update!(x, xdot, xdotdot, p, pdot, lambda, H_x, H_p, ϵ)
            interpolate_path!(x, alpha, s)
            _sgmam_refresh!(xdot, p, lambda, x, H_p)
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
        _sgmam_refresh!(xdot, p, lambda, x, H_p)
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
            _sgmam_refresh!(xdot, p, lambda, x, H_p)
            current_action = S_big
            stepsize = min(optimizer.stepsize_max, stepsize * optimizer.grow)
        else
            copyto!(x, x_start)
            _sgmam_refresh!(xdot, p, lambda, x, H_p)
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

function update!(x, xdot, xdotdot, p, pdot, lambda, H_x, H_p, ϵ)
    # xdot, p, and lambda are assumed to be pre-computed and consistent with x
    # (via _sgmam_refresh! or central_diff! + update_p!)
    central_diff!(pdot, p)
    Hx = H_x(x, p)

    # implicit update
    central_diff!(xdotdot, xdot)

    return update_x!(x, lambda, pdot, xdotdot, Hx, ϵ)
end

function update_x!(x, λ, p′, x′′, Hx, ϵ)
    # each dof has same lambda, but possibly different H_pp^{-1}
    Nx, Nt = size(x)
    xa = x[:, 1]
    xb = x[:, end]

    idxc = 2:(Nt - 1)

    # Tridiagonal system is shared across dof
    d = 1 .+ 2 .* ϵ .* λ[2:(end - 1)] .^ 2
    du = -ϵ .* λ[2:(end - 2)] .^ 2
    dl = -ϵ .* λ[3:(end - 1)] .^ 2
    T = LinearAlgebra.Tridiagonal(dl, d, du)

    # One cache per thread (solve! is not thread-safe on a shared cache)
    nthreads = Threads.nthreads()
    if hasmethod(Threads.nthreads, Tuple{Symbol})
        nthreads = Threads.nthreads(:default) + Threads.nthreads(:interactive)
    end
    caches = map(1:nthreads) do _
        b = zeros(eltype(x), length(d))
        prob = LinearProblem(T, b)
        init(
            prob,
            LUFactorization();
            alias = SciMLBase.LinearAliasSpecifier(; alias_A = true, alias_b = true),
        )
    end

    Threads.@threads for dof in 1:Nx
        cache = caches[Threads.threadid()]
        rhs = cache.b

        @inbounds for k in 1:length(rhs)
            i = k + 1
            rhs[k] = x[dof, i] + ϵ * (λ[i] * p′[dof, i] + Hx[dof, i] - λ[i]^2 * x′′[dof, i])
        end
        rhs[1] += ϵ * λ[2]^2 * xa[dof]
        rhs[end] += ϵ * λ[end - 1]^2 * xb[dof]

        solve!(cache)
        x[dof, idxc] .= cache.u
    end
    return nothing
end

function update_p!(p, lambda, x, xdot, H_p)
    # Alternative: Direct computation, only correct for quadratic Hamiltonian in p,
    # where H_pp does not depend on p
    b_ = H_p(x, 0 * x)
    lambda .= sqrt.(sum(b_ .^ 2; dims = 1) ./ sum(xdot .^ 2; dims = 1))
    lambda[1] = 0
    lambda[end] = 0
    p .= (lambda .* xdot .- b_)
    return nothing
end

function central_diff!(xdot, x)
    # ̇xₙ = 0.5(xₙ₊₁ - xₙ₋₁) central finite difference
    @views xdot[:, 2:(end - 1)] .= 0.5 .* (x[:, 3:end] .- x[:, 1:(end - 2)])
    return nothing
end

function _sgmam_refresh!(xdot, p, lambda, x, H_p)
    central_diff!(xdot, x)
    update_p!(p, lambda, x, xdot, H_p)
    return nothing
end

# Freidlin-Wentzell action of an instanton path discretized over a uniform arclength
# parameter `s ∈ [0, 1]`. With `xdot` produced by `central_diff!` (which does not
# divide by Δs), `xdot[:, i] ≈ Δs · dϕ/ds`, so `dot(xdot, p) ≈ ∫ p · dϕ`. On the
# zero-energy shell `H = p·f + |p|²/2 = 0`, the identity `p · dϕ = (|p|²/2) dt`
# holds along the instanton, so `∫ p · dϕ = ½ ∫ |p|² dt = S_FW`. No additional `/2`.
FW_action(xdot, p) = dot(xdot, p)
