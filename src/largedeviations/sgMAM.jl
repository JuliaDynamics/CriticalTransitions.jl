"""
    FreidlinWentzellHamiltonian{IIP, D, Hx, Hp, AF}

A discretized representation of the Freidlin-Wentzell Hamiltonian associated with an SDE
``dX = b(X)\\, dt + \\sqrt{\\varepsilon}\\, \\sigma(X)\\, dW``. The Hamiltonian is

```math
H(x, p) = \\langle b(x),\\, p \\rangle + \\tfrac{1}{2} \\langle p,\\, a(x)\\, p \\rangle,
```

where ``a(x) = \\sigma(x) \\sigma(x)^\\top``. The struct stores the partial derivatives
``H_x`` and ``H_p`` of the Hamiltonian and, for state-dependent noise, an optional
callable `a_func(x)` returning the diffusion tensor ``a(x)``. The iteration dispatches
on the return type of `a_func`:

* `a_func === nothing` → additive noise ``a \\equiv I``; the implicit step uses a single
  tridiagonal solve shared across all degrees of freedom (cheapest case).
* `a_func(x)::AbstractVector` → *diagonal* state-dependent ``a(x)``; the implicit step
  decouples and uses a tridiagonal solve **per** degree of freedom.
* `a_func(x)::AbstractMatrix` → *general* invertible ``a(x)``; the implicit step solves a
  sparse block-tridiagonal system of size ``(N-2)D \\times (N-2)D``.

Construct via:

* `FreidlinWentzellHamiltonian(ds)` for a `CoupledSDEs` / `CoupledODEs`. For
  multiplicative noise the constructor auto-selects the diagonal or general
  representation based on whether ``\\sigma(x)\\sigma(x)^\\top`` is diagonal at the
  initial state, and adds the state-dependent ``\\frac{1}{2} p^\\top (\\partial_x a) p``
  term to ``H_x`` via finite differences on ``a``.
* `FreidlinWentzellHamiltonian{IIP, D}(H_x, H_p; a_func = nothing)` for a hand-rolled
  Hamiltonian; the user-supplied `H_x` must already include the
  ``\\frac{1}{2} p^\\top (\\partial_x a) p`` term when `a_func` is state-dependent.

!!! note "Deprecated alias"
    `ExtendedPhaseSpace` is retained as a deprecated alias for backwards compatibility.
"""
struct FreidlinWentzellHamiltonian{IIP, D, Hx, Hp, AF}
    H_x::Hx
    H_p::Hp
    a_func::AF
    noise_shape::Symbol

    function FreidlinWentzellHamiltonian(ds::ContinuousTimeDynamicalSystem)
        if ds isa CoupledSDEs
            proper_MAM_system(ds)
        end

        f = dynamic_rule(ds)
        jac = jacobian(ds)
        param = current_parameters(ds)

        a_func, noise_shape = if ds isa CoupledSDEs
            _normalized_a_shape_func(ds)
        else
            nothing, :additive
        end

        function H_x(x, p)
            Hx = similar(x)
            for idx in 1:size(x, 2)
                xi = x[:, idx]
                pi = p[:, idx]
                jax = jac(xi, param, 0.0)
                @inbounds for idc in 1:size(x, 1)
                    Hx[idc, idx] = dot(jax[:, idc], pi)
                end
                if a_func !== nothing
                    Dn = size(x, 1)
                    h = max(sqrt(eps(eltype(xi))), 1.0e-8)
                    sample = a_func(xi)
                    is_diag = _is_diagonal_a(sample)
                    # Single preallocated unit-vector buffer reused across all j.
                    e = zeros(eltype(xi), Dn)
                    for j in 1:Dn
                        fill!(e, 0)
                        e[j] = h
                        ap = a_func(xi .+ e)
                        am = a_func(xi .- e)
                        if is_diag
                            # diagonal: ∂a_k/∂x_j; H_x_j += (1/2) Σ_k p_k² ∂a_k/∂x_j
                            dla = (_diag_view(ap) .- _diag_view(am)) ./ (2 * h)
                            Hx[j, idx] += 0.5 * dot(pi .^ 2, dla)
                        else
                            # general: ∂a_{kl}/∂x_j; H_x_j += (1/2) p^T (∂a/∂x_j) p
                            Hx[j, idx] += 0.5 * dot(pi, (ap .- am) ./ (2 * h), pi)
                        end
                    end
                end
            end
            return Hx
        end
        function H_p(x, p)
            Hp = similar(x)
            for idx in 1:size(x, 2)
                xi = x[:, idx]
                pi = p[:, idx]
                b = f(xi, param, 0.0)
                if a_func === nothing
                    Hp[:, idx] = pi .+ b
                else
                    a_i = a_func(xi)
                    if _is_diagonal_a(a_i)
                        Hp[:, idx] = b .+ _diag_view(a_i) .* pi
                    else
                        Hp[:, idx] = b .+ a_i * pi
                    end
                end
            end
            return Hp
        end
        return new{
            isinplace(ds), dimension(ds), typeof(H_x), typeof(H_p), typeof(a_func),
        }(H_x, H_p, a_func, noise_shape)
    end
    function FreidlinWentzellHamiltonian{IIP, D}(
            H_x::Function, H_p::Function; a_func = nothing,
        ) where {IIP, D}
        noise_shape = a_func === nothing ? :additive : :general
        return new{IIP, D, typeof(H_x), typeof(H_p), typeof(a_func)}(
            H_x, H_p, a_func, noise_shape,
        )
    end
end

Base.@deprecate_binding ExtendedPhaseSpace FreidlinWentzellHamiltonian

function prettyprint(mlp::FreidlinWentzellHamiltonian{IIP, D}) where {IIP, D}
    noise = mlp.noise_shape == :additive ? "additive" :
        mlp.noise_shape == :diagonal ? "diagonal multiplicative" :
        "general multiplicative"
    return "Freidlin-Wentzell Hamiltonian on $D-dimensional state space ($noise noise), containing $(IIP ? "in-place" : "out-of-place") H_x and H_p"
end

Base.show(io::IO, mlp::FreidlinWentzellHamiltonian) = print(io, prettyprint(mlp))

"""
$(TYPEDSIGNATURES)

Performs the simplified geometric Minimal Action Method (sgMAM) on the given system `sys`.

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
    H_p, H_x, a_func = sys.H_p, sys.H_x, sys.a_func

    Nx, Nt = size(x_initial)
    s = range(0; stop = 1, length = Nt)
    x, p, pdot, xdot, lambda, alpha = init_allocation(x_initial, Nt)
    xdotdot = zeros(size(xdot))
    p_zero = zero(xdot)

    x_prev = similar(x)

    # Ensure a consistent starting path for action comparisons
    interpolate_path!(x, alpha, s)
    _sgmam_refresh!(xdot, p, lambda, x, H_p, a_func, p_zero)
    initial_action = FW_action(xdot, p)

    function try_step!(ϵ)
        update!(x, xdot, xdotdot, p, pdot, lambda, H_x, H_p, ϵ, a_func)
        interpolate_path!(x, alpha, s)
        _sgmam_refresh!(xdot, p, lambda, x, H_p, a_func, p_zero)
        return FW_action(xdot, p)
    end
    save!() = copyto!(x_prev, x)
    function restore!()
        copyto!(x, x_prev)
        return _sgmam_refresh!(xdot, p, lambda, x, H_p, a_func, p_zero)
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
        optimizer::GeometricGradient = GeometricGradient(; stepsize = 1.0e3);
        kwargs...,
    )
    return minimize_simple_geometric_action(
        sys, Matrix(Matrix(x_initial)'), optimizer; kwargs...
    )
end
function minimize_simple_geometric_action(
        sys::ContinuousTimeDynamicalSystem,
        x_initial::Matrix{<:Real},
        optimizer::GeometricGradient = GeometricGradient(; stepsize = 1.0e3);
        kwargs...,
    )
    return minimize_simple_geometric_action(
        FreidlinWentzellHamiltonian(sys), Matrix(Matrix(x_initial)'), optimizer; kwargs...
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

function update!(x, xdot, xdotdot, p, pdot, lambda, H_x, H_p, ϵ, a_func = nothing)
    # xdot, p, and lambda are assumed to be pre-computed and consistent with x
    # (via _sgmam_refresh! or central_diff! + update_p!)
    central_diff!(pdot, p)
    Hx = H_x(x, p)

    # implicit update
    central_diff!(xdotdot, xdot)

    return update_x!(x, lambda, pdot, xdotdot, Hx, ϵ, a_func; path = x)
end

function update_x!(x, λ, p′, x′′, Hx, ϵ, a_func = nothing; path = x)
    # Three implicit-step variants, picked by the type of `a_func(x)`:
    # * a_func === nothing → additive a = I; single tridiagonal shared across DOFs.
    # * a_func(x) :: AbstractVector → diagonal a(x); per-DOF tridiagonal with coefficient
    #   λ_i² / a_kk(x_i). DOFs decouple in the implicit step.
    # * a_func(x) :: AbstractMatrix → general invertible a(x); the equation
    #   x_new[:, i] - ϵ λ_i² a^{-1}(x_i) (x_new[:, i+1] - 2 x_new[:, i] + x_new[:, i-1]) = R_i
    #   is multiplied by a(x_i) to avoid matrix inverses, giving a sparse block-tridiagonal
    #   system with diagonal blocks (a(x_i) + 2ϵλ_i² I) and off-diagonal blocks (-ϵλ_i² I).
    Nx, Nt = size(x)
    xa = x[:, 1]
    xb = x[:, end]
    idxc = 2:(Nt - 1)

    if a_func === nothing
        _update_x_additive!(x, λ, p′, x′′, Hx, ϵ, xa, xb, idxc, Nx, Nt)
    else
        sample = a_func(view(path, :, 2))
        if _is_diagonal_a(sample)
            _update_x_diag!(x, λ, p′, x′′, Hx, ϵ, xa, xb, idxc, Nx, Nt, a_func)
        else
            _update_x_general!(x, λ, p′, x′′, Hx, ϵ, xa, xb, idxc, Nx, Nt, a_func)
        end
    end
    return nothing
end

function _update_x_additive!(x, λ, p′, x′′, Hx, ϵ, xa, xb, idxc, Nx, Nt)
    d = 1 .+ 2 .* ϵ .* λ[2:(end - 1)] .^ 2
    du = -ϵ .* λ[2:(end - 2)] .^ 2
    dl = -ϵ .* λ[3:(end - 1)] .^ 2
    T = LinearAlgebra.Tridiagonal(dl, d, du)

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

function _update_x_diag!(x, λ, p′, x′′, Hx, ϵ, xa, xb, idxc, Nx, Nt, a_func)
    a_at = Matrix{eltype(x)}(undef, Nx, Nt)
    for i in 1:Nt
        a_at[:, i] .= _diag_view(a_func(x[:, i]))
    end

    Threads.@threads for dof in 1:Nx
        α = [ϵ * λ[i]^2 / a_at[dof, i] for i in 2:(Nt - 1)]
        d = 1 .+ 2 .* α
        du = -α[1:(end - 1)]
        dl = -α[2:end]
        T_dof = LinearAlgebra.Tridiagonal(dl, d, du)
        rhs = zeros(eltype(x), length(d))
        @inbounds for k in 1:length(rhs)
            i = k + 1
            rhs[k] = x[dof, i] + ϵ * (
                λ[i] * p′[dof, i] + Hx[dof, i] - (λ[i]^2 / a_at[dof, i]) * x′′[dof, i]
            )
        end
        rhs[1] += α[1] * xa[dof]
        rhs[end] += α[end] * xb[dof]
        x[dof, idxc] .= T_dof \ rhs
    end
    return nothing
end

function _update_x_general!(x, λ, p′, x′′, Hx, ϵ, xa, xb, idxc, Nx, Nt, a_func)
    N_in = Nt - 2
    n = N_in * Nx

    # Precompute a(x_i) at each interior grid point.
    A_blocks = Vector{Matrix{eltype(x)}}(undef, N_in)
    for i_in in 1:N_in
        A_blocks[i_in] = Matrix(a_func(x[:, i_in + 1]))
    end

    # Build sparse block-tridiagonal: diagonal blocks (a(x_i) + 2ϵλ_i² I),
    # off-diagonal blocks -ϵ λ_i² I.
    nnz_est = N_in * Nx * Nx + 2 * (N_in - 1) * Nx
    Iv = Vector{Int}(undef, 0); sizehint!(Iv, nnz_est)
    Jv = Vector{Int}(undef, 0); sizehint!(Jv, nnz_est)
    Vv = Vector{eltype(x)}(undef, 0); sizehint!(Vv, nnz_est)
    r = zeros(eltype(x), n)

    for i_in in 1:N_in
        i = i_in + 1
        base = (i_in - 1) * Nx
        A_i = A_blocks[i_in]
        λi2 = λ[i]^2
        # Diagonal block (a(x_i) + 2 ϵ λ_i² I)
        for j in 1:Nx, k in 1:Nx
            v = A_i[j, k] + (j == k ? 2 * ϵ * λi2 : zero(eltype(x)))
            push!(Iv, base + j); push!(Jv, base + k); push!(Vv, v)
        end
        # Off-diagonal block to the left: -ϵ λ_i² I
        if i_in > 1
            baseL = (i_in - 2) * Nx
            for j in 1:Nx
                push!(Iv, base + j); push!(Jv, baseL + j); push!(Vv, -ϵ * λi2)
            end
        end
        # Off-diagonal block to the right: -ϵ λ_i² I
        if i_in < N_in
            baseR = i_in * Nx
            for j in 1:Nx
                push!(Iv, base + j); push!(Jv, baseR + j); push!(Vv, -ϵ * λi2)
            end
        end
        # RHS: A_i (x[:, i] + ϵ(λ_i p'_i + Hx_i)) - ϵ λ_i² x''_i + boundary contributions
        rhs_block = A_i * (x[:, i] .+ ϵ .* (λ[i] .* p′[:, i] .+ Hx[:, i])) .-
                        ϵ * λi2 .* x′′[:, i]
        if i_in == 1
            rhs_block .+= ϵ * λi2 .* xa
        end
        if i_in == N_in
            rhs_block .+= ϵ * λi2 .* xb
        end
        @inbounds for j in 1:Nx
            r[base + j] = rhs_block[j]
        end
    end

    M = sparse(Iv, Jv, Vv, n, n)
    sol = M \ r
    for i_in in 1:N_in
        base = (i_in - 1) * Nx
        @inbounds for j in 1:Nx
            x[j, i_in + 1] = sol[base + j]
        end
    end
    return nothing
end

function update_p!(p, lambda, x, xdot, H_p, a_func = nothing, p_zero = zero(x))
    # Closed-form on the zero-energy shell H(x, p) = 0 (eqs. 40-41 of
    # Grafke-Schäfer-Vanden-Eijnden 2017):
    #   λ = |b|_a / |xdot|_a    with ⟨u, v⟩_a := ⟨u, a^{-1} v⟩
    #   p = a^{-1}(λ xdot - b)
    fill!(p_zero, 0)
    b_ = H_p(x, p_zero)
    Nt = size(x, 2)
    if a_func === nothing
        lambda .= sqrt.(sum(b_ .^ 2; dims = 1) ./ sum(xdot .^ 2; dims = 1))
        lambda[1] = 0
        lambda[end] = 0
        p .= (lambda .* xdot .- b_)
    else
        # Cache a(x_i) once per grid point to avoid calling `a_func` twice (once for
        # λ, once for p) — `a_func` is often the expensive part of the iteration.
        sample = a_func(view(x, :, 1))
        if _is_diagonal_a(sample)
            a_cache = Vector{Vector{eltype(x)}}(undef, Nt)
            @inbounds for i in 1:Nt
                a_cache[i] = collect(_diag_view(a_func(x[:, i])))
                a_i = a_cache[i]
                num = sum(b_[:, i] .^ 2 ./ a_i)
                den = sum(xdot[:, i] .^ 2 ./ a_i)
                lambda[1, i] = den > 0 ? sqrt(num / den) : zero(eltype(x))
            end
            lambda[1] = 0
            lambda[end] = 0
            @inbounds for i in 1:Nt
                a_i = a_cache[i]
                p[:, i] .= (lambda[1, i] .* xdot[:, i] .- b_[:, i]) ./ a_i
            end
        else
            a_cache = Vector{typeof(sample)}(undef, Nt)
            @inbounds for i in 1:Nt
                a_cache[i] = a_func(x[:, i])
                A_i = a_cache[i]
                A_inv_b = A_i \ b_[:, i]
                A_inv_xdot = A_i \ xdot[:, i]
                num = dot(b_[:, i], A_inv_b)
                den = dot(xdot[:, i], A_inv_xdot)
                lambda[1, i] = den > 0 ? sqrt(num / den) : zero(eltype(x))
            end
            lambda[1] = 0
            lambda[end] = 0
            @inbounds for i in 1:Nt
                A_i = a_cache[i]
                p[:, i] .= A_i \ (lambda[1, i] .* xdot[:, i] .- b_[:, i])
            end
        end
    end
    return nothing
end

function central_diff!(xdot, x)
    # ̇xₙ = 0.5(xₙ₊₁ - xₙ₋₁) central finite difference
    @views xdot[:, 2:(end - 1)] .= 0.5 .* (x[:, 3:end] .- x[:, 1:(end - 2)])
    return nothing
end

function _sgmam_refresh!(xdot, p, lambda, x, H_p, a_func = nothing, p_zero = zero(x))
    central_diff!(xdot, x)
    update_p!(p, lambda, x, xdot, H_p, a_func, p_zero)
    return nothing
end

"""
    FW_action(xdot, p)

Trapezoidal quadrature of ``\\int_0^1 \\langle p, \\dot\\varphi \\rangle\\, ds`` on the discretized path
used by the sgMAM iteration. On the zero-energy shell ``H(\\varphi, p) = 0`` this equals the
Freidlin-Wentzell action ``(1/2) \\int_0^T \\langle \\dot\\varphi - b,\\, a^{-1}(\\dot\\varphi - b) \\rangle\\, dt``.

`xdot` here is the *unnormalized* central difference `(x[:, i+1] - x[:, i-1]) / 2`
produced by [`central_diff!`](@ref), which already absorbs the arc-length measure `ds`, so no
further division is performed.
"""
FW_action(xdot, p) = dot(xdot, p)
