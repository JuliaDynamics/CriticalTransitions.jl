"""
$(TYPEDEF)

Multiple-shooting BVP optimizer for the Freidlin-Wentzell instanton on a
[`FreidlinWentzellHamiltonian`](@ref). Integrates the arclength-reparametrized Hamilton
equations
```math
\\frac{\\mathrm{d}\\varphi}{\\mathrm{d}s} = \\alpha\\,H_p,\\qquad
\\frac{\\mathrm{d}p}{\\mathrm{d}s}       = -\\alpha\\,H_x,\\qquad
\\alpha = L / \\|H_p\\|
```
on `s ∈ [0, 1]` with path length `L` a Newton unknown. Boundary states are parameterized by
the unstable / stable eigenvectors of the Hamiltonian Jacobian at each fixed-point endpoint.
The heteroclinic must not cross a drift fixed point in its interior; through-saddle
problems must be user-split into `attractor → saddle` legs and the actions summed.

# Fields
$(TYPEDFIELDS)
"""
struct MultipleShooting{S, J, T} <: GMAMOptimizer
    """Number of shooting segments (≥ 2)."""
    nshoots::Int
    """ODE integrator passed to the per-segment `solve`."""
    ode_solver::S
    """Override for the Newton solver; `nothing` builds a default sparse-AD `NewtonRaphson`."""
    nlsolve::J
    """Newton absolute tolerance."""
    abstol::T
    """Newton relative tolerance."""
    reltol::T
    """Anchor magnitude: distance of `y_0` from `(xa, 0)` in 2D-space."""
    eps_lin::T
    """Warning threshold on `max_s |H(φ(s), p(s))|` along the converged path."""
    invariant_tol::T
    """Maximum Newton iterations."""
    maxiters::Int
end

function MultipleShooting(;
        nshoots::Int = 10,
        ode_solver = OrdinaryDiffEqLowOrderRK.BS3(),
        nlsolve = nothing,
        abstol::Real = 1.0e-8,
        reltol::Real = 1.0e-6,
        eps_lin::Real = 1.0e-8,
        invariant_tol::Real = 1.0e-6,
        maxiters::Int = 50,
    )
    nshoots ≥ 2 || throw(ArgumentError("MultipleShooting requires nshoots ≥ 2; got $nshoots"))
    T = promote_type(
        typeof(float(abstol)), typeof(float(reltol)),
        typeof(float(eps_lin)), typeof(float(invariant_tol)),
    )
    return MultipleShooting{typeof(ode_solver), typeof(nlsolve), T}(
        nshoots, ode_solver, nlsolve,
        T(abstol), T(reltol), T(eps_lin), T(invariant_tol), maxiters,
    )
end

# The per-segment ODE must integrate ~3 digits tighter than the Newton tolerance so that
# integration error stays below the BVP residual floor that Newton is chasing.
const _SEG_TOL_FACTOR = 1.0e-3

function _drift(H::FreidlinWentzellHamiltonian{IIP, D}, x::AbstractVector) where {IIP, D}
    p0 = zeros(eltype(x), D)
    Hp = H.H_p(reshape(x, D, 1), reshape(p0, D, 1))
    return Hp[:, 1]
end

function _drift_jacobian(H::FreidlinWentzellHamiltonian{IIP, D}, x::AbstractVector) where {IIP, D}
    return ForwardDiff.jacobian(y -> _drift(H, y), x)
end

function _hamiltonian_value(H::FreidlinWentzellHamiltonian, x::AbstractVector, p::AbstractVector)
    b = _drift(H, x)
    a = collect(H.a(x))
    return LinearAlgebra.dot(b, p) + 0.5 * LinearAlgebra.dot(p, a * p)
end

function _assert_fixed_point(H::FreidlinWentzellHamiltonian, x::AbstractVector, side::Symbol)
    b = _drift(H, x)
    tol = sqrt(eps(real(eltype(x))))
    if LinearAlgebra.norm(b) > tol
        throw(ArgumentError(
            "MultipleShooting: $side endpoint x = $(collect(x)) is not a " *
            "fixed-point endpoint of the drift (‖b(x)‖ = $(LinearAlgebra.norm(b)) > √eps). " *
            "MultipleShooting requires both endpoints to be hyperbolic fixed points; " *
            "use GeometricGradient (gMAM/sgMAM) for free endpoints."
        ))
    end
    J = _drift_jacobian(H, x)
    ev = LinearAlgebra.eigvals(J)
    if any(abs(real(λ)) < tol for λ in ev)
        throw(ArgumentError(
            "MultipleShooting: $side endpoint x = $(collect(x)) is non-hyperbolic " *
            "(drift Jacobian has eigenvalue with |Re| ≈ 0; eigenvalues = $ev). " *
            "Either perturb the endpoint slightly off the bifurcation or use GeometricGradient."
        ))
    end
    return nothing
end

function _build_linearization(
        H::FreidlinWentzellHamiltonian{IIP, D}, x::AbstractVector, side::Symbol,
    ) where {IIP, D}
    T = float(real(eltype(x)))
    J = T.(_drift_jacobian(H, x))
    A = collect(T, H.a(collect(x)))
    M = zeros(T, 2D, 2D)
    @views begin
        M[1:D, 1:D] .= J
        M[1:D, (D + 1):(2D)] .= A
        M[(D + 1):(2D), (D + 1):(2D)] .= -J'
    end
    Sf = LinearAlgebra.schur(M)
    select = if side === :outgoing
        [real(λ) > 0 for λ in Sf.values]
    elseif side === :incoming
        [real(λ) < 0 for λ in Sf.values]
    else
        throw(ArgumentError("side must be :outgoing or :incoming, got $side"))
    end
    LinearAlgebra.ordschur!(Sf, select)
    U = Matrix{T}(Sf.Z[:, 1:D])
    Λ = Matrix{T}(Sf.T[1:D, 1:D])
    xstar_aug = vcat(collect(T, x), zeros(T, D))
    return (xstar_aug = xstar_aug, U = U, Λ = Λ)
end

struct MultipleShootingWorkspace{IIP, D, HT, T, ODE, NLS, LINA, LINB, SP}
    H::HT
    nshoots::Int
    grid::Vector{T}
    ode_solver::ODE
    nlsolve::NLS
    lin_a::LINA
    lin_b::LINB
    sparsity::SP
    abstol::T
    reltol::T
    eps_lin::T
    invariant_tol::T
    maxiters::Int
end

_ws_D(::MultipleShootingWorkspace{IIP, D}) where {IIP, D} = D

function _arclength_rhs!(dy, y, p_params, s)
    H = p_params.H
    D = p_params.D
    L = p_params.L
    φ_mat = reshape(view(y, 1:D), D, 1)
    p_mat = reshape(view(y, (D + 1):(2D)), D, 1)
    Hp = H.H_p(φ_mat, p_mat)
    Hx = H.H_x(φ_mat, p_mat)
    norm_Hp_sq = zero(eltype(y))
    @inbounds for k in 1:D
        norm_Hp_sq += Hp[k, 1]^2
    end
    norm_Hp = sqrt(norm_Hp_sq)
    α = L / max(norm_Hp, sqrt(eps(real(eltype(y)))))
    @inbounds for k in 1:D
        dy[k] = α * Hp[k, 1]
        dy[D + k] = -α * Hx[k, 1]
    end
    return nothing
end

function _integrate_segment(ws::MultipleShootingWorkspace{IIP, D}, y_in, s_a, s_b, L) where {IIP, D}
    params = (H = ws.H, D = D, L = L)
    prob = SciMLBase.ODEProblem(_arclength_rhs!, collect(y_in), (s_a, s_b), params)
    sol = SciMLBase.solve(
        prob, ws.ode_solver;
        abstol = ws.abstol * _SEG_TOL_FACTOR,
        reltol = ws.reltol * _SEG_TOL_FACTOR,
        save_everystep = false, save_start = false, dense = false,
    )
    return sol.u[end]
end

function _unpack_unknowns(z, D, nseg)
    c_a = view(z, 1:D)
    offset = D
    int_len = 2D * (nseg - 1)
    interior_flat = view(z, (offset + 1):(offset + int_len))
    offset += int_len
    c_b = view(z, (offset + 1):(offset + D))
    L = z[end]
    return c_a, interior_flat, c_b, L
end

_residual_size(D::Int, nseg::Int) = 2D * nseg + 1

# Return a view of node i's (φ, p) state for i ∈ 0:nseg. Endpoints come from the
# linearization parameterization; interior nodes are direct slices of `interior_flat`.
function _node_state(D, nseg, y0, yend, interior_flat, i)
    if i == 0
        return y0
    elseif i == nseg
        return yend
    else
        return view(interior_flat, ((i - 1) * 2D + 1):(i * 2D))
    end
end

function _shooting_residual!(F, z, ws::MultipleShootingWorkspace{IIP, D}) where {IIP, D}
    nseg = ws.nshoots
    c_a, interior_flat, c_b, L = _unpack_unknowns(z, D, nseg)

    diff_lin_a = ws.lin_a.U * c_a
    y0 = ws.lin_a.xstar_aug .+ diff_lin_a
    yend = ws.lin_b.xstar_aug .+ ws.lin_b.U * c_b

    fidx = 0
    for i in 1:nseg
        y_in = _node_state(D, nseg, y0, yend, interior_flat, i - 1)
        y_target = _node_state(D, nseg, y0, yend, interior_flat, i)
        y_end = _integrate_segment(ws, y_in, ws.grid[i], ws.grid[i + 1], L)
        @inbounds for k in 1:(2D)
            F[fidx + k] = y_end[k] - y_target[k]
        end
        fidx += 2D
    end

    F[fidx + 1] = sum(abs2, diff_lin_a) - ws.eps_lin^2
    return nothing
end

function _build_jac_sparsity(D::Int, nseg::Int)
    nz = _residual_size(D, nseg)
    rows = Int[]
    cols = Int[]
    ca_cols = 1:D
    interior_col(i) = (D + (i - 1) * 2D + 1):(D + i * 2D)
    cb_cols = (D + 2D * (nseg - 1) + 1):(D + 2D * (nseg - 1) + D)
    L_col = nz

    function push_block!(rrange, crange)
        for r in rrange, c in crange
            push!(rows, r); push!(cols, c)
        end
    end

    for i in 1:nseg
        rrange = ((i - 1) * 2D + 1):(i * 2D)
        if i == 1
            push_block!(rrange, ca_cols)
        else
            push_block!(rrange, interior_col(i - 1))
        end
        if i == nseg
            push_block!(rrange, cb_cols)
        else
            push_block!(rrange, interior_col(i))
        end
        push_block!(rrange, L_col:L_col)
    end
    push_block!((2D * nseg + 1):(2D * nseg + 1), ca_cols)

    return SparseArrays.sparse(rows, cols, ones(Float64, length(rows)), nz, nz)
end

function _build_workspace(
        H::FreidlinWentzellHamiltonian{IIP, D},
        x_init::AbstractMatrix,
        opt::MultipleShooting,
    ) where {IIP, D}
    size(x_init, 1) == D || throw(ArgumentError(
        "x_init has $(size(x_init, 1)) rows; expected D = $D for this FreidlinWentzellHamiltonian."
    ))
    xa = x_init[:, 1]
    xb = x_init[:, end]
    _assert_fixed_point(H, xa, :outgoing)
    _assert_fixed_point(H, xb, :incoming)
    lin_a = _build_linearization(H, xa, :outgoing)
    lin_b = _build_linearization(H, xb, :incoming)
    T = promote_type(float(eltype(x_init)), typeof(opt.abstol))
    nseg = opt.nshoots
    grid = collect(range(zero(T), one(T); length = nseg + 1))
    sparsity = _build_jac_sparsity(D, nseg)
    return MultipleShootingWorkspace{
        IIP, D,
        typeof(H), T,
        typeof(opt.ode_solver), typeof(opt.nlsolve),
        typeof(lin_a), typeof(lin_b),
        typeof(sparsity),
    }(
        H, nseg, grid, opt.ode_solver, opt.nlsolve, lin_a, lin_b,
        sparsity,
        T(opt.abstol), T(opt.reltol),
        T(opt.eps_lin), T(opt.invariant_tol), opt.maxiters,
    )
end

# Seed (c, L)-style unstable/stable mode coefficient from the configuration tangent at the
# endpoint, projected onto the linearization subspace and rescaled so the 2D-norm equals
# `eps_lin` (matches the anchor residual at the outgoing side; arbitrary at the incoming
# side, but a consistent magnitude is a decent Newton starting point).
function _project_endpoint(lin, x_near, eps_lin::T) where {T}
    xstar = lin.xstar_aug[1:length(x_near)]
    tangent = vcat(collect(T, x_near .- xstar), zero(x_near))
    c_raw = lin.U' * tangent
    norm_lin = LinearAlgebra.norm(lin.U * c_raw)
    scale = norm_lin > eps(T) ? eps_lin / norm_lin : eps_lin
    return T.(c_raw .* scale)
end

function _initial_guess_unknowns(ws::MultipleShootingWorkspace{IIP, D}, x_init) where {IIP, D}
    T = eltype(ws.grid)
    nseg = ws.nshoots
    N = size(x_init, 2)
    c_a = _project_endpoint(ws.lin_a, T.(x_init[:, min(2, N)]), ws.eps_lin)
    c_b = _project_endpoint(ws.lin_b, T.(x_init[:, max(N - 1, 1)]), ws.eps_lin)
    interior = zeros(T, 2D * (nseg - 1))
    for i in 1:(nseg - 1)
        idx = clamp(round(Int, (i / nseg) * (N - 1)) + 1, 1, N)
        φ_i = T.(x_init[:, idx])
        b_i = _drift(ws.H, φ_i)
        a_i = collect(T, ws.H.a(φ_i))
        p_i = a_i \ (-2 .* b_i)
        @inbounds for k in 1:D
            interior[(i - 1) * 2D + k] = φ_i[k]
            interior[(i - 1) * 2D + D + k] = p_i[k]
        end
    end
    L0 = zero(T)
    @inbounds for i in 1:(N - 1)
        L0 += LinearAlgebra.norm(T.(x_init[:, i + 1]) .- T.(x_init[:, i]))
    end
    return vcat(c_a, interior, c_b, [T(max(L0, ws.eps_lin))])
end

function _solve_shooting(ws::MultipleShootingWorkspace{IIP, D}, z0) where {IIP, D}
    nz = _residual_size(D, ws.nshoots)
    f! = (F, z, _) -> _shooting_residual!(F, z, ws)
    jac_proto = SparseArrays.SparseMatrixCSC{eltype(z0), Int}(ws.sparsity)
    coloring = SparseMatrixColorings.coloring(
        ws.sparsity,
        SparseMatrixColorings.ColoringProblem{:nonsymmetric, :column}(),
        GreedyColoringAlgorithm(),
    )
    colorvec = SparseMatrixColorings.column_colors(coloring)
    nf = SciMLBase.NonlinearFunction(
        f!;
        resid_prototype = zeros(eltype(z0), nz),
        jac_prototype = jac_proto,
        colorvec = colorvec,
    )
    prob = SciMLBase.NonlinearProblem(nf, z0, nothing)
    alg = ws.nlsolve === nothing ?
        NonlinearSolveFirstOrder.NewtonRaphson(
            autodiff = AutoForwardDiff(),
            linsolve = LinearSolve.UMFPACKFactorization(),
        ) :
        ws.nlsolve
    sol = SciMLBase.solve(prob, alg; abstol = ws.abstol, reltol = ws.reltol, maxiters = ws.maxiters)
    return sol.u, LinearAlgebra.norm(sol.resid), sol.stats.nsteps, sol.retcode
end

function _sample_path(ws::MultipleShootingWorkspace{IIP, D}, z, N::Int) where {IIP, D}
    T = eltype(ws.grid)
    nseg = ws.nshoots
    c_a, interior_flat, c_b, L = _unpack_unknowns(z, D, nseg)
    y0 = ws.lin_a.xstar_aug .+ ws.lin_a.U * c_a
    yend = ws.lin_b.xstar_aug .+ ws.lin_b.U * c_b
    path = zeros(T, D, N)
    pmat = zeros(T, D, N)
    arclength = collect(range(zero(T), one(T); length = N))
    H_inv_max = zero(T)
    for j in 1:N
        s = arclength[j]
        i = clamp(searchsortedfirst(ws.grid, s) - 1, 1, nseg)
        s_start = ws.grid[i]
        y_in = _node_state(D, nseg, y0, yend, interior_flat, i - 1)
        y_at_s = s == s_start ? collect(y_in) : _integrate_segment(ws, y_in, s_start, s, L)
        @inbounds for k in 1:D
            path[k, j] = y_at_s[k]
            pmat[k, j] = y_at_s[D + k]
        end
        Hval = _hamiltonian_value(ws.H, view(y_at_s, 1:D), view(y_at_s, (D + 1):(2D)))
        H_inv_max = max(H_inv_max, abs(Hval))
    end
    return path, pmat, arclength, T(L), H_inv_max
end

function minimize_geometric_action(
        sys::FreidlinWentzellHamiltonian,
        x_init::AbstractMatrix,
        opt::MultipleShooting;
        show_progress::Bool = false,
    )
    N = size(x_init, 2)
    ws = _build_workspace(sys, x_init, opt)
    z0 = _initial_guess_unknowns(ws, x_init)
    z_sol, resnorm, niter, retcode = _solve_shooting(ws, z0)
    if retcode != SciMLBase.ReturnCode.Success
        throw(ErrorException(
            "MultipleShooting: NonlinearSolve did not converge (retcode = $retcode, " *
            "residual = $resnorm, iterations = $niter). Try increasing nshoots, " *
            "warm-starting from a GeometricGradient solve, or tightening abstol."
        ))
    end
    path_mat, p_mat, arclength, L, H_inv_max = _sample_path(ws, z_sol, N)
    if H_inv_max > ws.invariant_tol
        @warn "MultipleShooting: H=0 invariant violated" H_invariant_max=H_inv_max threshold=ws.invariant_tol
    end
    b_fn = x -> _drift(sys, x)
    A_at = x -> inv(collect(sys.a(x)))
    v_buf = similar(path_mat)
    integrand_buf = zeros(eltype(path_mat), N)
    action = _geometric_action_from_drift!(b_fn, path_mat, one(eltype(path_mat)), A_at, v_buf, integrand_buf)
    show_progress && @info "MultipleShooting converged" residual=resnorm iterations=niter path_length=L action=action
    return MinimumActionPath(
        StateSpaceSet(Matrix(path_mat')), action;
        λ = L, generalized_momentum = Matrix(p_mat'),
    )
end

function minimize_geometric_action(
        ::CoupledSDEs, ::Matrix, ::MultipleShooting; kwargs...,
    )
    throw(ArgumentError(
        "MultipleShooting acts on FreidlinWentzellHamiltonian, not CoupledSDEs. " *
        "The shooting method operates on the LDT-derived deterministic Hamiltonian system, " *
        "not the original stochastic system. Convert explicitly:\n\n" *
        "    H = FreidlinWentzellHamiltonian(sys)\n" *
        "    minimize_geometric_action(H, x_init, MultipleShooting(...))\n"
    ))
end
