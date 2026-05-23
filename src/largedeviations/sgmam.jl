"""
    minimize_geometric_action(sys::FreidlinWentzellHamiltonian, x_initial, optimizer = GeometricGradient(; stepsize = 1e3); kwargs...)

Runs the simplified geometric Minimum Action Method (sgMAM, [grafke_long_2017](@cite)) on
the [`FreidlinWentzellHamiltonian`](@ref) `sys`, returning the instanton (minimizer of the
Freidlin-Wentzell action) that connects the two endpoints of `x_initial`.

sgMAM works directly in the doubled ``(x, p)`` phase space and minimizes the symplectic
action ``\\int p \\cdot \\mathrm{d}x`` under the zero-energy constraint ``H(x, p) = 0``.
Here ``p = a(x)^{-1}(\\dot x - b(x))`` is the conjugate momentum, i.e. the (dimensionless)
**noise force** that the Wiener process must supply to drive the trajectory over the barrier.

## Arguments
* `sys::FreidlinWentzellHamiltonian`: the Hamiltonian, typically obtained from a
  `CoupledSDEs` via `FreidlinWentzellHamiltonian(ds)`.
* `x_initial::Matrix{T}` (or `StateSpaceSet`): initial guess for the instanton, shape
  `D × Nt` with state points in columns. Endpoints `x_initial[:, 1]` and `x_initial[:, end]`
  are held fixed; the interior is reparameterized to uniform arclength on the first iteration.
* `optimizer`: step-size control, either
  - [`GeometricGradient`](@ref) (default; backtracking projected gradient), or
  - [`AdaptiveGeometricGradient`](@ref) (multi-phase probe variant; more robust on
    underdamped systems at higher per-iteration cost).

## Keyword arguments
* `maxiters::Int = 1000`: maximum outer iterations.
* `abstol::Real = NaN`, `reltol::Real = NaN`: convergence on the absolute / relative
  action change between accepted iterations. `NaN` disables that criterion.
* `verbose::Bool = false`: log convergence and backtracking events.
* `show_progress::Bool = false`: show a `ProgressMeter` bar.

Returns a [`MinimumActionPath`](@ref)
"""
function minimize_geometric_action(
        sys::FreidlinWentzellHamiltonian,
        x_initial::Matrix{T},
        optimizer::GeometricGradient = GeometricGradient(; stepsize = 1.0e3);
        maxiters::Int = 1000, show_progress::Bool = false, verbose::Bool = false,
        abstol::Real = NaN, reltol::Real = NaN,
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
        maxiters::Int = 1000, show_progress::Bool = false, verbose::Bool = false,
        abstol::Real = NaN, reltol::Real = NaN,
    ) where {T}
    Nt = size(x_initial, 2)
    cache = build_sgmam_cache(sys, x_initial, Nt)
    return _minimize_geometric_action_inner_adaptive(
        sys, x_initial, optimizer, cache;
        maxiters, show_progress, verbose, abstol, reltol,
    )
end

# Shared initialization: allocate buffers, arclength-reparameterize, refresh λ/p/xdot.
function _init_sgmam_state(sys, x_initial, cache)
    Nt = size(x_initial, 2)
    s = range(0; stop = 1, length = Nt)
    x, p, pdot, xdot, lambda, alpha = init_allocation(x_initial, Nt)
    xdotdot = zeros(size(xdot))
    interpolate_path!(x, alpha, s)
    _sgmam_refresh!(xdot, p, lambda, x, sys, cache)
    return x, p, pdot, xdot, xdotdot, lambda, alpha, s
end

function _minimize_geometric_action_inner(
        sys, x_initial, optimizer::GeometricGradient, cache;
        maxiters, show_progress, verbose, abstol, reltol,
    )
    x, p, pdot, xdot, xdotdot, lambda, alpha, s = _init_sgmam_state(sys, x_initial, cache)
    x_prev = similar(x)
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
        sys, x_initial, optimizer::AdaptiveGeometricGradient, cache;
        maxiters, show_progress, verbose, abstol, reltol,
    )
    x, p, pdot, xdot, xdotdot, lambda, alpha, s = _init_sgmam_state(sys, x_initial, cache)
    x_start = similar(x)
    x_big_result = similar(x)
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
            stepsize <= optimizer.stepsize_min && break
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
            break
        end
        next!(
            progress; step = n,
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
