"""
$(TYPEDSIGNATURES)

Performs the simplified geometric Minimal Action Method (sgMAM) on the given system `sys`.
Our implementation is only valid for additive noise.

This method computes the optimal path in the phase space of a Hamiltonian system that
minimizes the Freidlin-Wentzell action. The Hamiltonian functions `H_x` and `H_p` define
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
