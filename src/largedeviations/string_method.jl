"""
    $(TYPEDSIGNATURES)

Compute the string method for a given system using the method of [e_string_2007](@citet).

The string method is an iterative algorithm used to find minimum energy path (MEP) between two points in phase space. It works by evolving a discretized path (string) according to the system's drift while maintaining equal arc-length parametrization between points.

This implementation allows for computation between arbitrary points, not just stable fixed points.

## References

- [e_string_2007](@cite)

# Arguments
- `sys::ExtendedPhaseSpace`: The doubled phase space system for which the string method is computed
- `x_initial`: Initial path discretized as a matrix where each column represents a point on the path
- `ϵ::Real`: Step size for the evolution step
- `alg`: SciML integrator algorithm (e.g. `Euler()`, `Tsit5()`). Defaults to `Euler()`.
- `iterations::Int64`: Maximum number of iterations for path convergence
- `show_progress::Bool`: Whether to display a progress meter during computation

# Returns
- `x`: The final converged path representing the MEP
"""
function string_method(
    sys::Union{ExtendedPhaseSpace,Function},
    x_initial::Matrix;
    ϵ::Real=1e-1,
    alg=Euler(),
    iterations::Int64=1000,
    show_progress::Bool=false,
)
    Nx, Nt = size(x_initial)
    s = range(0; stop=1, length=Nt)
    x, alpha = init_allocation_string(x_initial, Nt)

    integ = _init_string_integrator(sys, x; ϵ, alg)

    progress = Progress(iterations; dt=0.5, enabled=show_progress)
    for i in 1:iterations
        _sciml_step!(x, integ, ϵ)
        # reset initial and final points to allow for string computation
        # between points that are not stable fixed points
        x[:, 1] = x_initial[:, 1]
        x[:, end] = x_initial[:, end]

        # reparameterize to arclength
        interpolate_path!(x, alpha, s)

        next!(progress; showvalues=[("iterations", i)])
    end
    return StateSpaceSet(x')
end
"""
    $(TYPEDSIGNATURES)

Compute the string method for a given system using the method of [e_string_2007](@citet).

The string method is an iterative algorithm used to find minimum energy path (MEP) between two points in phase space. It works by evolving a discretized path (string) according to the system's drift while maintaining equal arc-length parametrization between points.

This implementation allows for computation between arbitrary points, not just stable fixed points.

## References

- [e_string_2007](@cite)

# Arguments
- `sys::CoupledSDEs`: The system for which the string method is computed
- `x_initial`: Initial path discretized as a matrix where each column represents a point on the path
- `ϵ::Real`: Step size for the evolution step
- `alg`: SciML integrator algorithm (e.g. `Euler()`, `Tsit5()`). Defaults to `Euler()`.
- `iterations::Int64`: Maximum number of iterations for path convergence
- `show_progress::Bool`: Whether to display a progress meter during computation

# Returns
- `x`: The final converged path representing the MEP
"""
function string_method(sys::ContinuousTimeDynamicalSystem, init; kwargs...)
    b(x) = drift(sys, x)
    return string_method(b, init; kwargs...)
end

function string_method(
    b::Union{ExtendedPhaseSpace,Function},
    x_initial::StateSpaceSet{D};
    ϵ::Real=1e-1,
    alg=Euler(),
    iterations::Int64=1000,
    show_progress::Bool=false,
) where {D}
    x_initial_m, Nt = _string_matrix(x_initial)
    s = range(0; stop=1, length=Nt)

    x, alpha = init_allocation_string(x_initial_m, Nt)
    integ = _init_string_integrator(b, x; ϵ, alg)

    progress = Progress(iterations; dt=0.5, enabled=show_progress)
    for i in 1:iterations
        _sciml_step!(x, integ, ϵ)

        # reset initial and final points to allow for string computation
        # between points that are not stable fixed points
        x[:, 1] = x_initial_m[:, 1]
        x[:, end] = x_initial_m[:, end]

        # reparameterize to arclength
        interpolate_path!(x, alpha, s)

        next!(progress; showvalues=[("iterations", i)])
    end
    return StateSpaceSet(x')
end

function _string_matrix(x::StateSpaceSet{D}) where {D}
    m = Matrix(x)
    Nt = length(x)
    if size(m) == (Nt, D)
        return permutedims(m), Nt
    elseif size(m) == (D, Nt)
        return m, Nt
    else
        throw(
            ArgumentError(
                "Unexpected Matrix(StateSpaceSet) size $(size(m)); expected $(Nt)×$(D) or $(D)×$(Nt).",
            ),
        )
    end
end

function _string_matrix(y::AbstractMatrix, D::Int, Nt::Int)
    if size(y) == (D, Nt)
        return y
    elseif size(y) == (Nt, D)
        return permutedims(y)
    else
        throw(
            ArgumentError(
                "Unexpected matrix size $(size(y)); expected $(D)×$(Nt) or $(Nt)×$(D)."
            ),
        )
    end
end

function _string_matrix(y::StateSpaceSet{D}, ::Val{D}, Nt::Int) where {D}
    m = Matrix(y)
    return _string_matrix(m, D, Nt)
end

function _hp_mode(sys::ExtendedPhaseSpace, x::Matrix)
    D, Nt = size(x)

    y_m = sys.H_p(x, zeros(size(x)))
    if y_m isa AbstractMatrix && size(y_m) == size(x)
        return :matrix
    end

    x_sss = StateSpaceSet(x')
    p0_sss = StateSpaceSet(zeros(Nt, D))
    y_s = sys.H_p(x_sss, p0_sss)
    if y_s isa StateSpaceSet
        y_sm = Matrix(y_s)
        if size(y_sm) == (Nt, D) || size(y_sm) == (D, Nt)
            return :sss
        end
    end

    throw(
        ArgumentError(
            "`ExtendedPhaseSpace.H_p` must return either a $(D)×$(Nt) Matrix when given a Matrix state, or a $(Nt)×$(D) StateSpaceSet when given a StateSpaceSet state.",
        ),
    )
end

function _init_string_integrator(sys::Function, x::Matrix; ϵ::Real, alg)
    _, Nt = size(x)
    function f!(du, u, p, t)
        fill!(du, 0)
        @inbounds @views for j in 2:(Nt - 1)
            du[:, j] .= sys(u[:, j])
        end
        return nothing
    end

    prob = SciMLBase.ODEProblem{true}(f!, x, (0.0, float(ϵ)))
    return SciMLBase.init(
        prob,
        alg;
        adaptive=false,
        dt=float(ϵ),
        save_everystep=false,
        save_start=false,
        save_end=false,
        dense=false,
    )
end

function _init_string_integrator(sys::ExtendedPhaseSpace, x::Matrix; ϵ::Real, alg)
    D, Nt = size(x)
    mode = _hp_mode(sys, x)

    p0_m = zeros(size(x))
    x_sss = StateSpaceSet(x')
    p0_sss = StateSpaceSet(zeros(Nt, D))

    function f!(du, u, p, t)
        if mode === :matrix
            du .= sys.H_p(u, p0_m)
        else
            x_sss = StateSpaceSet(u')
            du .= _string_matrix(sys.H_p(x_sss, p0_sss), Val(D), Nt)
        end

        @views begin
            du[:, 1] .= 0
            du[:, end] .= 0
        end
        return nothing
    end

    prob = SciMLBase.ODEProblem{true}(f!, x, (0.0, float(ϵ)))
    return SciMLBase.init(
        prob,
        alg;
        adaptive=false,
        dt=float(ϵ),
        save_everystep=false,
        save_start=false,
        save_end=false,
        dense=false,
    )
end

function _sciml_step!(x::Matrix, integ, ϵ::Real)
    SciMLBase.reinit!(integ, x; t0=0.0, tf=float(ϵ), erase_sol=true, reset_dt=false)
    SciMLBase.step!(integ, float(ϵ), true)
    copyto!(x, integ.u)
    return nothing
end

function init_allocation_string(x_initial, Nt)
    # preallocate
    x = deepcopy(x_initial)
    alpha = zeros(Nt)
    return x, alpha
end
