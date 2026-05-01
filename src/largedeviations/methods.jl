abstract type GMAMOptimizer end

"""
$(TYPEDEF)

Optimizer configuration for the (s)gMAM projected-gradient update with built-in
backtracking step-size control.

By default, backtracking is **enabled** (`max_backtracks=10`): each iteration tries the
current step size and, if the action increases, shrinks it by `shrink` and retries up to
`max_backtracks` times. On accepted steps the step size grows by `grow` for the next
iteration. This makes the algorithm insensitive to the initial step size choice: small
starting values are grown automatically, and large values are safely reduced by backtracking.
Set `max_backtracks=0` to disable backtracking and use a fixed step size.

# Fields
$(TYPEDFIELDS)

# Keyword constructors
$(METHODLIST)
"""
struct GeometricGradient{T<:Real} <: GMAMOptimizer
    """Initial step size for the projected gradient update."""
    stepsize::T
    """Step-size shrink factor on rejected steps (backtracking)."""
    shrink::T
    """Step-size growth factor after accepted steps (backtracking)."""
    grow::T
    """Maximum number of backtracking attempts per iteration."""
    max_backtracks::Int
    """Lower clamp for step size."""
    stepsize_min::T
    """Upper clamp for step size."""
    stepsize_max::T
end

function GeometricGradient(;
    stepsize::Real=1e-1,
    shrink::Real=0.5,
    grow::Real=1.1,
    max_backtracks::Int=10,
    stepsize_min::Real=1e-12,
    stepsize_max::Real=Inf,
)
    T = promote_type(
        typeof(float(stepsize)),
        typeof(float(shrink)),
        typeof(float(grow)),
        typeof(float(stepsize_min)),
        typeof(float(stepsize_max)),
    )
    return GeometricGradient{T}(
        T(stepsize), T(shrink), T(grow), max_backtracks, T(stepsize_min), T(stepsize_max)
    )
end

"""
    backtracking_optimize!(optimizer, try_step!, save!, restore!, initial_action; kwargs...)

Generic backtracking iteration loop shared by sgMAM and gMAM.

## Arguments
  - `optimizer::GeometricGradient`: step-size control parameters
  - `try_step!(ϵ) -> action`: perform one update at step size `ϵ` and return the resulting action
  - `save!()`: save current state (called before each backtracking attempt)
  - `restore!()`: revert to saved state (called on rejected steps and retries)
  - `initial_action`: action value at the starting path

## Keyword arguments
  - `maxiters`, `abstol`, `reltol`, `verbose`, `show_progress`: same as the solver functions

Returns `(final_action, final_stepsize)`.
"""
function backtracking_optimize!(
    optimizer::GeometricGradient,
    try_step!,
    save!,
    restore!,
    initial_action::Real;
    maxiters::Int=1000,
    abstol::Real=NaN,
    reltol::Real=NaN,
    verbose::Bool=false,
    show_progress::Bool=false,
)
    backtracking = optimizer.max_backtracks > 0
    ntries = backtracking ? optimizer.max_backtracks + 1 : 1
    stepsize = optimizer.stepsize
    current_action = oftype(stepsize, initial_action)
    consecutive_failures = 0
    max_consecutive_failures = 5

    progress = Progress(maxiters; dt=0.5, enabled=show_progress)
    for i in 1:maxiters
        S_old = current_action
        ϵ_try = if backtracking
            clamp(stepsize, optimizer.stepsize_min, optimizer.stepsize_max)
        else
            stepsize
        end
        accepted = !backtracking

        backtracking && save!()
        for try_idx in 1:ntries
            if backtracking && try_idx > 1
                restore!()
            end

            S_trial = oftype(stepsize, try_step!(ϵ_try))

            if !backtracking || (isfinite(S_trial) && S_trial <= S_old)
                accepted = true
                current_action = S_trial
                if backtracking
                    stepsize = min(optimizer.stepsize_max, ϵ_try * optimizer.grow)
                end
                break
            end

            if backtracking
                ϵ_try *= optimizer.shrink
                if ϵ_try < optimizer.stepsize_min
                    break
                end
            end
        end

        if backtracking
            if accepted
                consecutive_failures = 0
            else
                restore!()
                current_action = S_old
                stepsize = max(optimizer.stepsize_min, ϵ_try)
                consecutive_failures += 1
                if stepsize <= optimizer.stepsize_min
                    verbose &&
                        @info "Step-size search stalled at stepsize_min=$(optimizer.stepsize_min) after $i iterations."
                    break
                end
                if consecutive_failures >= max_consecutive_failures
                    verbose &&
                        @info "Step-size search stalled: $consecutive_failures consecutive backtrack failures after $i iterations."
                    break
                end
            end
        end

        abs_change = accepted ? abs(current_action - S_old) : Inf
        rel_change = if accepted
            current_action == 0 ? abs_change : abs_change / abs(current_action)
        else
            Inf
        end

        if accepted && (
            (isfinite(abstol) && abs_change < abstol) ||
            (isfinite(reltol) && rel_change < reltol)
        )
            verbose &&
                @info "Converged after $i iterations with abs=$abs_change, rel=$rel_change"
            break
        end
        next!(
            progress;
            showvalues=[
                ("iterations", i),
                ("Stol", round(rel_change; sigdigits=3)),
                ("stepsize", round(stepsize; sigdigits=3)),
            ],
        )
    end
    return current_action, stepsize
end
