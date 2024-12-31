"""
A structure representing a system with Hamiltonian functions H_x and H_p.

This system operates in an extended phase space where the Hamiltonian is assumed to be
quadratic in the extended momentum. The phase space coordinates `x` are doubled to
form a 2n-dimensional extended phase space.
"""
struct SgmamSystem
    H_x::Function
    H_p::Function
end

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
function sgmam(
    sys::SgmamSystem,
    x_initial::Matrix{<:Real};
    ϵ::Float64=1e-1,
    iterations::Int64=1000,
    show_progress::Bool=false,
    reltol::Float64=NaN,
)
    H_p, H_x = sys.H_p, sys.H_x

    Nx, Nt = size(x_initial)
    s = range(0; stop=1, length=Nt)
    x, p, pdot, xdot, lambda, alpha = init_allocation(x_initial, Nt)

    S = CircularBuffer{Float64}(2)
    fill!(S, Inf)

    progress = Progress(iterations; dt=0.5, enabled=show_progress)
    for i in 1:iterations
        update!(x, xdot, p, pdot, lambda, H_x, H_p, ϵ)

        # reparametrize to arclength
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
function sgmam(sys, x_initial::StateSpaceSet; kwargs...)
    return sgmam(sys, Matrix(Matrix(x_initial)'); kwargs...)
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
