"""
A structure representing a system with Hamiltonian functions H_x and H_p.

This system operates in an extended phase space where the Hamiltonian is assumed to be
quadratic in the extended momentum. The phase space coordinates `x` are doubled to
form a 2n-dimensional extended phase space.
"""
struct ExtendedHamiltonianSystem{IIP,D,Hx,Hp}
    H_x::Hx
    H_p::Hp

    function ExtendedHamiltonianSystem(ds::ContinuousTimeDynamicalSystem)
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
    function ExtendedHamiltonianSystem{IIP,D}(H_x::Function, H_p::Function) where {IIP,D}
        return new{IIP,D,typeof(H_x),typeof(H_p)}(H_x, H_p)
    end
end

function prettyprint(mlp::ExtendedHamiltonianSystem{IIP,D}) where {IIP,D}
    return "Doubled $D-dimensional phase space containing $(IIP ? "in-place" : "out-of-place") functions"
end

Base.show(io::IO, mlp::ExtendedHamiltonianSystem) = print(io, prettyprint(mlp))

"""
$(TYPEDSIGNATURES)

Performs the simplified geometric Minimal Action Method (sgMAM) on the given system `sys`.
Our implementation is only valid for additive noise.

This method computes the optimal path in the phase space of a Hamiltonian system that
minimizes the Freidlin–Wentzell action. The Hamiltonian functions `H_x` and `H_p` define
the system's dynamics in a doubled phase. The initial state `x_initial` is evolved
iteratively using constrained gradient descent with step size parameter `ϵ` over a specified
number of iterations. The method can display a progress meter and will stop early if the
relative tolerance `reltol` is achieved.

The function returns a tuple containing the final state, the action value,
the Lagrange multipliers, the momentum, and the state derivatives. The implementation is
based on the work of [Grafke et al. (2019)](https://homepages.warwick.ac.uk/staff/T.Grafke/simplified-geometric-minimum-action-method-for-the-computation-of-instantons.html.
).
"""
function simple_geometric_min_action_method(
    sys::ExtendedHamiltonianSystem,
    x_initial::Matrix{T};
    ϵ::Real=1e-1,
    iterations::Int64=1000,
    show_progress::Bool=false,
    reltol::Real=NaN,
) where {T}
    H_p, H_x = sys.H_p, sys.H_x

    Nx, Nt = size(x_initial)
    s = range(0; stop=1, length=Nt)
    x, p, pdot, xdot, lambda, alpha = init_allocation(x_initial, Nt)

    S = CircularBuffer{T}(2)
    fill!(S, Inf)

    progress = Progress(iterations; dt=0.5, enabled=show_progress)
    for i in 1:iterations
        update!(x, xdot, p, pdot, lambda, H_x, H_p, ϵ)

        # reparameterize to arclength
        interpolate_path!(x, alpha, s)
        push!(S, FW_action(xdot, p))

        tol = abs(S[end] - S[1]) / S[end]
        if tol < reltol
            @info "Converged after $i iterations with $tol"
            break
        end
        next!(progress; showvalues=[("iterations", i), ("Stol", round(tol; sigdigits=3))])
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
        ExtendedHamiltonianSystem(sys), Matrix(Matrix(x_initial)'); kwargs...
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

function update!(x, xdot, p, pdot, lambda, H_x, H_p, ϵ)
    central_diff!(xdot, x)

    update_p!(p, lambda, x, xdot, H_p)

    central_diff!(pdot, p)
    Hx = H_x(x, p)

    # implicit update
    xdotdot = zeros(size(xdot))
    central_diff!(xdotdot, xdot)

    return update_x!(x, lambda, pdot, xdotdot, Hx, ϵ)
end

function update_x!(x, λ, p′, x′′, Hx, ϵ)
    # each dof has same lambda, but possibly different H_pp^{-1}
    Nx, Nt = size(x)
    xa = x[:, 1]
    xb = x[:, end]

    # rhs = zeros(Nt)
    idxc = 2:(Nt - 1)
    Threads.@threads for dof in 1:Nx
        rhs = @. (
            x[dof, idxc] +
            ϵ * (λ[idxc] * p′[dof, idxc] + Hx[dof, idxc] - λ[idxc]^2 * x′′[dof, idxc])
        )
        rhs[1] += ϵ * λ[2]^2 * xa[dof]
        rhs[end] += ϵ * λ[end - 1]^2 * xb[dof]

        A = spdiagm( # spdiagm makes it significantly faster
            0 => 1 .+ 2 .* ϵ .* λ[2:(end - 1)] .^ 2,
            1 => -ϵ .* λ[2:(end - 2)] .^ 2,
            -1 => -ϵ .* λ[3:(end - 1)] .^ 2,
        )
        prob = LinearProblem(A, rhs)
        x[dof, 2:(end - 1)] .= solve(prob, KLUFactorization()).u
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
