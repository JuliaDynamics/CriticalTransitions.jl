"""
$(TYPEDSIGNATURES)

Performs the simplified geometric Minimal Action Method (sgMAM) on the given Hamiltonian
system `sys`.

The method computes the optimal instanton path that minimizes the Freidlin-Wentzell action
in the doubled `(x, p)` phase space via constrained gradient descent. Returns a
[`MinimumActionPath`](@ref) containing the final path, action value, Lagrange multipliers
(`.λ`), generalized momentum (`.generalized_momentum`), and path velocity (`.path_velocity`).

The optional positional argument `optimizer` controls step-size adaptation. It defaults to
`GeometricGradient(; stepsize = 1e3)`. Pass `GeometricGradient(; max_backtracks = 0)` to use
a fixed step size, or `AdaptiveGeometricGradient(...)` for the multi-phase variant.

## Keyword arguments
  - `maxiters::Int = 1000`
  - `show_progress::Bool = false`
  - `verbose::Bool = false`
  - `abstol::Real = NaN`, `reltol::Real = NaN`: action-change convergence criteria
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
    Nt = size(x_initial, 2)
    cache = build_sgmam_cache(sys, x_initial, Nt)
    return _minimize_geometric_action_inner(
        sys, x_initial, optimizer, cache;
        maxiters, show_progress, verbose, abstol, reltol,
    )
end

function minimize_geometric_action(
        sys::FreidlinWentzellHamiltonian,
        x_initial::StateSpaceSet,
        optimizer::GMAMOptimizer = GeometricGradient(; stepsize = 1.0e3);
        kwargs...,
    )
    return minimize_geometric_action(sys, Matrix(Matrix(x_initial)'), optimizer; kwargs...)
end

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
    Nt = size(x_initial, 2)
    cache = build_sgmam_cache(sys, x_initial, Nt)
    return _minimize_geometric_action_inner_adaptive(
        sys, x_initial, optimizer, cache;
        maxiters, show_progress, verbose, abstol, reltol,
    )
end

function _minimize_geometric_action_inner(
        sys::FreidlinWentzellHamiltonian, x_initial::Matrix{T}, optimizer::GeometricGradient, cache;
        maxiters, show_progress, verbose, abstol, reltol,
    ) where {T}
    Nt = size(x_initial, 2)
    s = range(0; stop = 1, length = Nt)
    x, p, pdot, xdot, lambda, alpha = init_allocation(x_initial, Nt)
    xdotdot = zeros(size(xdot))
    x_prev = similar(x)

    interpolate_path!(x, alpha, s)
    _sgmam_refresh!(xdot, p, lambda, x, sys, cache)
    initial_action = FW_action(xdot, p)

    function try_step!(ϵ)
        update!(x, xdot, xdotdot, p, pdot, lambda, sys, ϵ, cache)
        interpolate_path!(x, alpha, s)
        _sgmam_refresh!(xdot, p, lambda, x, sys, cache)
        return FW_action(xdot, p)
    end
    save!() = copyto!(x_prev, x)
    function restore!()
        copyto!(x, x_prev)
        return _sgmam_refresh!(xdot, p, lambda, x, sys, cache)
    end

    current_action, _ = backtracking_optimize!(
        optimizer, try_step!, save!, restore!, initial_action;
        maxiters, abstol, reltol, verbose, show_progress,
    )
    return MinimumActionPath(
        StateSpaceSet(x'), current_action;
        λ = lambda, generalized_momentum = p, path_velocity = xdot,
    )
end

function _minimize_geometric_action_inner_adaptive(
        sys::FreidlinWentzellHamiltonian, x_initial::Matrix{T},
        optimizer::AdaptiveGeometricGradient, cache;
        maxiters, show_progress, verbose, abstol, reltol,
    ) where {T}
    Nt = size(x_initial, 2)
    s = range(0; stop = 1, length = Nt)
    x, p, pdot, xdot, lambda, alpha = init_allocation(x_initial, Nt)
    xdotdot = zeros(size(xdot))

    x_start = similar(x)
    x_big_result = similar(x)

    interpolate_path!(x, alpha, s)
    _sgmam_refresh!(xdot, p, lambda, x, sys, cache)
    Tϵ = typeof(optimizer.stepsize)
    current_action = Tϵ(FW_action(xdot, p))

    function run_probe!(ϵ, n)
        S = oftype(ϵ, NaN)
        for _ in 1:n
            update!(x, xdot, xdotdot, p, pdot, lambda, sys, ϵ, cache)
            interpolate_path!(x, alpha, s)
            _sgmam_refresh!(xdot, p, lambda, x, sys, cache)
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
        _sgmam_refresh!(xdot, p, lambda, x, sys, cache)
        ϵ_small = clamp(
            stepsize * optimizer.shrink, optimizer.stepsize_min, optimizer.stepsize_max
        )
        S_small = run_probe!(ϵ_small, n)

        big_ok = isfinite(S_big) && S_big <= S_prev
        small_ok = isfinite(S_small) && S_small <= S_prev
        accepted = true
        if small_ok && (!big_ok || S_small < S_big)
            current_action = S_small
            stepsize = max(optimizer.stepsize_min, stepsize * optimizer.shrink)
        elseif big_ok
            copyto!(x, x_big_result)
            _sgmam_refresh!(xdot, p, lambda, x, sys, cache)
            current_action = S_big
            stepsize = min(optimizer.stepsize_max, stepsize * optimizer.grow)
        else
            copyto!(x, x_start)
            _sgmam_refresh!(xdot, p, lambda, x, sys, cache)
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
        StateSpaceSet(x'), current_action;
        λ = lambda, generalized_momentum = p, path_velocity = xdot,
    )
end

function init_allocation(x_initial, Nt)
    x = deepcopy(x_initial)
    p = zeros(size(x))
    pdot = zeros(size(x))
    xdot = zeros(size(x))
    lambda = zeros(1, Nt)
    alpha = zeros(Nt)
    return x, p, pdot, xdot, lambda, alpha
end
