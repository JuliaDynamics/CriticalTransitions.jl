"""
A structure representing a extanded phase space system where your dissipative vector field is embedded in a doubled dimensional phase space. Given old phase space coordinates `x` of a vector field `f(x)`, we can define the canonical momenta `p`, such that the new phase space coordinates are `(x, p)`. The dynamics in this extended phase space are then governed by the Hamtiltonian system:

``H = p^2 + x \\dot f(x)``

Hence, this system operates in an extended phase space where the Hamiltonian is assumed to be quadratic in the extended momentum.

The struct `ExtendedPhaseSpace` holds the Hamilton's functions `H_x` and `H_p`.
"""
struct ExtendedPhaseSpace{IIP,D,Hx,Hp}
    H_x::Hx
    H_p::Hp

    function ExtendedPhaseSpace(ds::ContinuousTimeDynamicalSystem)
        if ds isa CoupledSDEs
            proper_sgMAM_system(ds)
        end

        f = dynamic_rule(ds)
        jac = jacobian(ds)
        param=current_parameters(ds)

        function H_x(x, p) # ℜ² → ℜ²
            Hx = similar(x)

            for idx in 1:size(x, 2)
                jax = jac(x[:, idx], param, 0.0)
                for idc in 1:size(x, 1)
                    Hx[idc, idx] = dot(jax[:, idc], p[:, idx])
                end
            end
            return Hx
        end
        function H_p(x, p) # ℜ² → ℜ²
            Hp = similar(x)

            for idx in 1:size(x, 2)
                Hp[:, idx] = p[:, idx] + f(x[:, idx], param, 0.0)
            end
            return Hp
        end
        return new{isinplace(ds),dimension(ds),typeof(H_x),typeof(H_p)}(H_x, H_p)
    end
    function ExtendedPhaseSpace{IIP,D}(H_x::Function, H_p::Function) where {IIP,D}
        return new{IIP,D,typeof(H_x),typeof(H_p)}(H_x, H_p)
    end
end

function prettyprint(mlp::ExtendedPhaseSpace{IIP,D}) where {IIP,D}
    return "Doubled $D-dimensional phase space containing $(IIP ? "in-place" : "out-of-place") functions"
end

Base.show(io::IO, mlp::ExtendedPhaseSpace) = print(io, prettyprint(mlp))

"""
$(TYPEDSIGNATURES)

Performs the simplified geometric Minimal Action Method (sgMAM) on the given system `sys`.
Our implementation is only valid for additive noise.

This method computes the optimal path in the phase space of a Hamiltonian system that
minimizes the Freidlin–Wentzell action. The Hamiltonian functions `H_x` and `H_p` define
the system's dynamics in a doubled phase. The initial state `x_initial` is evolved
iteratively using constrained gradient descent with step size parameter `stepsize` over a specified
number of iterations. The method can display a progress meter and will stop early if the
absolute tolerance `abstol` or relative tolerance `reltol` is achieved.

The function returns a [`MinimumActionPath`](@ref) containing the final path, the action value,
the Lagrange multipliers (`.λ`), the momentum (`.generalized_momentum`), and the state derivatives (`.path_velocity`).
The implementation is based on the work of [Grafke et al. (2019)](https://homepages.warwick.ac.uk/staff/T.Grafke/simplified-geometric-minimum-action-method-for-the-computation-of-instantons.html).

## Keyword arguments

  - `stepsize::Real=1e-1`: step size for gradient descent. Default: `0.1`
  - `maxiters::Int=1000`: maximum number of iterations before the algorithm stops
  - `show_progress::Bool=false`: if true, display a progress bar
  - `verbose::Bool=false`: if true, print additional output
  - `abstol::Real=NaN`: absolute tolerance for early stopping based on action change
  - `reltol::Real=NaN`: relative tolerance for early stopping based on action change
"""
function simple_geometric_min_action_method(
    sys::ExtendedPhaseSpace,
    x_initial::Matrix{T};
    stepsize::Real=1e-1,
    maxiters::Int=1000,
    show_progress::Bool=false,
    verbose::Bool=false,
    abstol::Real=NaN,
    reltol::Real=NaN,
) where {T}
    H_p, H_x = sys.H_p, sys.H_x

    Nx, Nt = size(x_initial)
    s = range(0; stop=1, length=Nt)
    x, p, pdot, xdot, lambda, alpha = init_allocation(x_initial, Nt)
    xdotdot = zeros(size(xdot))

    S = CircularBuffer{T}(2)
    fill!(S, Inf)

    progress = Progress(maxiters; dt=0.5, enabled=show_progress)
    for i in 1:maxiters
        update!(x, xdot, xdotdot, p, pdot, lambda, H_x, H_p, stepsize)

        # reparameterize to arclength
        interpolate_path!(x, alpha, s)
        push!(S, FW_action(xdot, p))

        abs_change = abs(S[end] - S[1])
        rel_change = S[end] == 0 ? abs_change : abs_change / abs(S[end])
        if (isfinite(abstol) && abs_change < abstol) ||
            (isfinite(reltol) && rel_change < reltol)
            verbose && @info "Converged after $i iterations with abs=$abs_change, rel=$rel_change"
            break
        end
        next!(
            progress;
            showvalues=[("iterations", i), ("Stol", round(rel_change; sigdigits=3))],
        )
    end
    return MinimumActionPath(
        StateSpaceSet(x'), S[end]; λ=lambda, generalized_momentum=p, path_velocity=xdot
    )
end
function simple_geometric_min_action_method(sys, x_initial::StateSpaceSet; kwargs...)
    return simple_geometric_min_action_method(sys, Matrix(Matrix(x_initial)'); kwargs...)
end
function simple_geometric_min_action_method(
    sys::ContinuousTimeDynamicalSystem, x_initial::Matrix{<:Real}; kwargs...
)
    return simple_geometric_min_action_method(
        ExtendedPhaseSpace(sys), Matrix(Matrix(x_initial)'); kwargs...
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
    central_diff!(xdot, x)

    update_p!(p, lambda, x, xdot, H_p)

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
            alias=SciMLBase.LinearAliasSpecifier(; alias_A=true, alias_b=true),
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
    lambda .= sqrt.(sum(b_ .^ 2; dims=1) ./ sum(xdot .^ 2; dims=1))
    lambda[1] = 0
    lambda[end] = 0
    p .= (lambda .* xdot .- b_)
    return nothing
end

function central_diff!(xdot, x)
    # ̇xₙ = 0.5(xₙ₊₁ - xₙ₋₁) central finite difference
    xdot[:, 2:(end - 1)] = 0.5 * (x[:, 3:end] - x[:, 1:(end - 2)])
    return nothing
end

FW_action(xdot, p) = sum(sum(xdot .* p; dims=1)) / 2

function proper_sgMAM_system(ds::CoupledSDEs)
    proper_MAM_system(ds)
    Σ = covariance_matrix(ds)
    return LinearAlgebra.isdiag(Σ) || throw(
        ArgumentError(
            "Simple geometric action is only defined for diagonal noise. The noise covariance matrix is not diagonal.",
        ),
    )
end
