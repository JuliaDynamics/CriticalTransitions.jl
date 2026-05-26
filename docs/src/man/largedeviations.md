# Large deviation theory

This section applies results of large deviation theory (LDT), particularly action minimization problems for computing most probable transition paths in stochastic dynamical systems driven by weak noise. For a description of the theory, see [freidlin_random_1998](@citet) and [borner_climate_2025](@cite). An overview of numerical methods applying LDT is given in [grafke_numerical_2019](@citet).

!!! info
    The methods in this section apply to ``D``-dimensional stochastic dynamical systems of the form
    ```math
    \text{d} \mathbf{x} = \mathbf{b} (\mathbf{x}) \text{d}t + \sigma \mathbf{\Sigma} \text{d}\mathbf{W}_t \,,
    ```
    where the drift field ``\mathbf{b}`` may be non-gradient but the noise term must consist of Gaussian noise (``\mathbf{W}_t`` is a ``D``-dimensional vector of independent standard Wiener processes) and a constant covariance matrix ``\mathbf{Q} = \mathbf{\Sigma}\mathbf{\Sigma}^\top``.

    This is a special case of the broader class of noise types supported by [`CoupledSDEs`](@ref).

!!! note "Noise convention used by the action functionals"
    **Starting point.** A [`CoupledSDEs`](@ref) exposes its noise via the diffusion
    function ``g`` and, for additive invertible noise, the noise rate covariance
    ```math
    \mathbf{\Sigma} \;=\; \texttt{covariance\_matrix(sys)} \;=\; g\,g^\top.
    ```
    For an SDE built as ``\mathrm{d}\mathbf{x} = \mathbf{b}\,\mathrm{d}t + \sigma\,\mathbf{\Sigma}_0\,\mathrm{d}\mathbf{W}_t``
    with ``\mathbf{\Sigma}_0\mathbf{\Sigma}_0^\top = \mathbf{Q}``, this gives ``\mathbf{\Sigma} = \sigma^2\,\mathbf{Q}``.
    The Freidlin–Wentzell rate function only depends on the *product* ``\sigma^2\mathbf{Q}``;
    the individual ``(\sigma, \mathbf{Q})`` split is not determined by the SDE alone
    (``(\sigma, \mathbf{Q}) \leftrightarrow (c\sigma, \mathbf{Q}/c^2)`` is the same physical system).

    **Convention.** The package picks a unique ``(\sigma, \mathbf{Q})`` pair from ``\mathbf{\Sigma}``
    by adopting the **trace convention** ``\mathrm{tr}(\mathbf{Q}) = D``:
    ```math
    \sigma_{\mathrm{eff}}^2 \;=\; \frac{\mathrm{tr}(\mathbf{\Sigma})}{D},
    \qquad
    \mathbf{Q}_{\mathrm{can}} \;=\; \frac{D}{\mathrm{tr}(\mathbf{\Sigma})}\,\mathbf{\Sigma}.
    ```
    The accessor [`noise_strength`](@ref)`(sys)` returns ``\sigma_{\mathrm{eff}}``. The trace is
    invariant under orthogonal changes of basis and reduces to the per-direction average
    noise variance, making this the natural choice among rotation-invariant normalizations.

    **Action functionals use the canonical covariance.** All action functionals in this
    section evaluate the FW integrand in the ``\mathbf{Q}_{\mathrm{can}}^{-1}`` metric:
    ```math
    \texttt{fw\_action(sys, path, time)} \;=\; \Phi_{FW}^{\mathbf{Q}_{\mathrm{can}}}[\phi] \;=\;
        \tfrac{1}{2}\!\int (\dot\phi - \mathbf{b})^\top \mathbf{Q}_{\mathrm{can}}^{-1}(\dot\phi - \mathbf{b})\,\mathrm{d}t.
    ```
    Two consequences:

    - The returned action is **independent of the `noise_strength` keyword** chosen at
      construction (both ``\sigma`` and the magnitude of ``\mathbf{Q}`` are absorbed into
      ``\sigma_{\mathrm{eff}}``).
    - The returned action is **invariant under orthogonal changes of basis**, matching
      the coordinate-independence of the FW action as a geometric quantity.

    **Converting to the user's parameterization.** If you built your system as
    ``\mathrm{d}\mathbf{x} = \mathbf{b}\,\mathrm{d}t + \sigma_{\mathrm{user}}\,\mathbf{\Sigma}_0\,\mathrm{d}\mathbf{W}_t``
    with ``\mathbf{Q}_{\mathrm{user}} = \mathbf{\Sigma}_0\mathbf{\Sigma}_0^\top`` and want the FW action in
    *your* metric ``\mathbf{Q}_{\mathrm{user}}^{-1}``, multiply the returned action by
    ``\mathrm{tr}(\mathbf{Q}_{\mathrm{user}})/D``:
    ```math
    \Phi_{FW}^{\mathbf{Q}_{\mathrm{user}}}[\phi] \;=\; \frac{\mathrm{tr}(\mathbf{Q}_{\mathrm{user}})}{D}\cdot
        \texttt{fw\_action(sys, path, time)}.
    ```
    The factor is ``1`` whenever ``\mathbf{Q}_{\mathrm{user}}`` is isotropic (``c\,\mathbf{I}`` for any
    ``c>0``) or already trace-normalized (``\mathrm{tr}(\mathbf{Q}_{\mathrm{user}}) = D``), in which
    case `fw_action` returns ``\Phi_{FW}^{\mathbf{Q}_{\mathrm{user}}}`` directly. All examples and
    tests in the package fall in this case.

    **`fw_action` vs `om_action`.** [`fw_action`](@ref) and [`geometric_action`](@ref) are
    fully ``\sigma``-independent under this convention. [`om_action`](@ref) takes
    `noise_strength` as an explicit positional argument because the Onsager–Machlup
    correction ``(\sigma^2/2)\int \nabla\cdot \mathbf{b}\,\mathrm{d}t`` is a *finite-noise*
    correction that scales with ``\sigma^2``; only the FW part of OM is ``\sigma``-independent.
    Pass the ``\sigma`` at which you want ``-\log P[\phi]`` evaluated (typically the same
    `noise_strength` you used at construction). As ``\sigma \to 0`` the OM correction
    vanishes and OM ``\to`` FW, recovering the leading-order LDP rate function.

## Action minimizers
Several methods have been proposed to calculate transition paths that minimize a given [action functional](@ref "Action functionals"). In the weak-noise limit, this minimum action path (or instanton) corresponds to the most probable transition path. While the minimum action method (MAM) is the most basic version, it is often beneficial to minimize the [geometric action](@ref "Geometric Freidlin-Wentzell action") via a time-independent version called gMAM. The problem can also be cast in a Hamiltonian form, implemented as simple gMAM (sgMAM), which can have numerical advantages.

These methods apply to non-gradient systems driven by Gaussian noise. In gradient systems, minimum action paths between attractors coincide with heteroclinic orbits, which can be computed via the so-called string method.

!!! info "Action minimization as an optimal control problem"
    All of the action minimizers listed below (MAM, gMAM, sgMAM) can equivalently be cast as **optimal control problems**: find a path ``\mathbf{x}`` and a control ``\mathbf{u}`` that minimize a Freidlin-Wentzell-type action
    ```math
    \int \| \mathbf{u}(t) - \mathbf{b}(\mathbf{x}(t)) \|^2 \, \text{d}t
    ```
    subject to the path dynamics ``\dot{\mathbf{x}}(t) = \mathbf{u}(t)`` and the boundary conditions ``\mathbf{x}(0) = \mathbf{x}_A``, ``\mathbf{x}(T) = \mathbf{x}_B`` (with ``T`` either fixed, as in MAM, or eliminated by reparameterization, as in gMAM/sgMAM). In this form the problem can be solved with dedicated optimal-control packages such as [`OptimalControl.jl`](https://control-toolbox.org/OptimalControl.jl); see the [control-toolbox.org MAM tutorial](https://control-toolbox.org/Tutorials.jl/stable/tutorial-mam.html) for a worked example on the Maier-Stein system.

To summarize, the following methods are available:

- Minimum action method [(MAM)](@ref "Minimum action method (MAM)")
- Geometric minimum action method [(gMAM)](@ref "Geometric minimum action method (gMAM)")
- Simple geometric minimum action method [(sgMAM)](@ref "Simple geometric minimum action method (sgMAM)")
- [String method](@ref)

All three action-minimizers (MAM, gMAM, sgMAM) require only **autonomous** noise and reject **rank-deficient** diffusion; they all support additive, diagonal-multiplicative, and general-matrix multiplicative diffusion. The Onsager-Machlup functional (selectable in MAM via `functional = "OM"`) is the one exception: it is implemented only for additive noise and throws otherwise. Pick by what you want out of the result:

| Method | Use when |
|---|---|
| **MAM** | The transition time $T$ is fixed and you want the action at that $T$ (or you specifically need the Onsager-Machlup functional, which only has a time-parameterized form). |
| **gMAM / sgMAM** | $T$ is unknown and you want a time-reparameterization-invariant instanton. Prefer **sgMAM** (Hamiltonian picture) when an analytic `FreidlinWentzellHamiltonian(H_x, H_p)` is available or when AD-based `jacobian(ds)` evaluations are cheap; prefer **gMAM** for the closer-to-textbook geometric-action formulation. |
| **String method** | You want the deterministic heteroclinic orbit (typical use: **gradient** systems, where it coincides with the instanton). In non-gradient systems it generally differs from the most probable transition path. |

#### Variants and extensions
The literature contains a number of extensions of MAM-type methods that may be relevant depending on the model class and numerical difficulties. These variants are not currently implemented in `CriticalTransitions.jl`, but serve as useful pointers:

- **tMAM / optimal linear time scaling**: avoids explicit optimization over the transition time by introducing an optimal linear time scaling; can be combined with adaptivity in time discretization [wan_tmam_2015](@citet).
- **Adaptive MAM**: uses a moving-mesh strategy to concentrate grid points in dynamically important portions of the path, improving efficiency and robustness [zhou_adaptive_mam_2008](@citet).
- **Non-Gaussian (jump / Lévy) noise**: for systems driven by jump noise, the rate function and path optimization problem differ from the Freidlin--Wentzell diffusive setting; see e.g. an optimal-control-based approach in [wei_most_likely_jumps_2023](@citet).
- **Multiplicative / state-dependent noise**: state-dependent diagonal and full-matrix multiplicative noise is supported by both `gMAM` (via `minimize_geometric_action(::CoupledSDEs, ...)`) and `sgMAM` (via `FreidlinWentzellHamiltonian(::CoupledSDEs)`). The diffusion tensor `a(x)` is classified once when the sgMAM cache is built (constant vs state-dependent, diagonal vs coupled) and the resulting inner loop dispatches at compile time on the cache type; see [grafke_small_random_2017](@citet) for the underlying Hamiltonian formulation. Rank-deficient (degenerate) noise is rejected at cache build; see [grafke_small_random_2017](@citet) for the rank-deficient extension that would be needed to support it.

### Minimum action method (MAM)
Minimization of the specified action functional using the optimization algorithm of `Optimization.jl`. See also [e_minimum_2004](@citet).

```@docs
minimize_action
```

### Geometric minimum action method (gMAM)
Minimization of the geometric action following [heymann_pathways_2008, heymann_geometric_2008](@citet). gMAM reformulates MAM to avoid double optimization of both the action and the transition time. It achieves this by using a [geometric action](@ref "Geometric Freidlin-Wentzell action") functional that is independent of the time parametrization of the path. This reparameterization invariance makes the method more robust and computationally efficient, particularly for multiscale systems.

```@docs
minimize_geometric_action
```

### Simple geometric minimum action method (sgMAM)
Simplified minimization of the geometric action following [grafke_long_2017](@citet).
The simple gMAM reduces the complexity of the original gMAM by requiring only first-order derivatives of the underlying Hamiltonian optimization formulation. This simplifies the numerical treatment and computational complexity.

The implementation below performs a constrained gradient descent on the Hamiltonian system, supporting autonomous diffusions with additive, diagonal-multiplicative, or general-matrix multiplicative noise.

#### Hamiltonian picture

The Freidlin-Wentzell rate function admits an equivalent Hamiltonian formulation with conjugate momentum ``p``:

```math
H(\varphi, p) = \langle b(\varphi), p \rangle + \tfrac{1}{2}\, \langle p,\, a(\varphi)\, p \rangle,
```

where ``a(x) = \sigma(x)\sigma(x)^{\top}``. The action along an instanton (where ``H \equiv 0``) reduces to ``S[\varphi] = \int_0^T \langle p, \mathrm{d}\varphi\rangle``, the form used by `minimize_geometric_action` when the input is a [`FreidlinWentzellHamiltonian`](@ref).

Compile-time dispatch on the shape of ``a(x)`` is driven by the type of `sys.a` (a `Base.Returns` wrapper marks constant-in-`x` diffusion) together with the type of `a(x_ref)` (a `LinearAlgebra.Diagonal` marks diagonal coupling). Classification happens once when the per-path cache is built; the resulting cache type then selects the diagonal or coupled inner loops in `update_p!`, `update_x!`, and `geometric_gradient_step!`.

Rank-deficient ``a(x)`` is rejected at cache build (`a` is probed at a reference state `x_ref` (the first path point) and at the ``2D`` neighbors `x_ref ± h·eₗ`).

```@docs
FreidlinWentzellHamiltonian
```


#### Performance notes (sgMAM)

sgMAM repeatedly evaluates `H_p(x, p)` and `H_x(x, p)` along a discretized path. If this allocates, sgMAM will be slow.

Key reason for performance differences:

- `FreidlinWentzellHamiltonian(ds::CoupledSDEs)` typically relies on `jacobian(ds)` (often automatic differentiation unless you provide an analytic Jacobian) and evaluates it pointwise along the path.
- A hardcoded `FreidlinWentzellHamiltonian(H_x, H_p)` with analytic expressions operating on the full `D×Nt` path matrix usually allocates far less.

Benchmark pattern:

```julia
using CriticalTransitions
using BenchmarkTools

sys_fast = FreidlinWentzellHamiltonian{false,2}(H_x, H_p)  # hardcoded analytic H_x/H_p

ds = CoupledSDEs(KPO, zeros(2), ())
sys_generic = FreidlinWentzellHamiltonian(ds)              # uses jacobian(ds)

opt = GeometricGradient(; stepsize=0.5)
@btime minimize_geometric_action($sys_fast,    $x_initial, $opt; maxiters=100, show_progress=false)
@btime minimize_geometric_action($sys_generic, $x_initial, $opt; maxiters=100, show_progress=false)
```

Aside: the same “vectorized + allocation-free inner loop” principle also tends to make [`string_method`](@ref) faster when used with `FreidlinWentzellHamiltonian`.

### Multiple shooting

The [`MultipleShooting`](@ref) optimizer treats the instanton as a boundary value problem on arclength-reparametrized Hamilton equations
```math
\frac{\mathrm{d}\varphi}{\mathrm{d}s} = \alpha\,H_p(\varphi, p),\qquad
\frac{\mathrm{d}p}{\mathrm{d}s} = -\alpha\,H_x(\varphi, p),\qquad
\alpha = L\,/\,\|H_p\|
```
on `s ∈ [0, 1]`, with the total path length `L` a Newton unknown so that `‖dφ/ds‖ ≡ L` pointwise. The boundary states `y_0, y_{N_seg}` are parameterized by the unstable / stable eigenvectors of the Hamiltonian Jacobian `M = [J A; 0 -Jᵀ]` at each fixed-point endpoint. The BVP is segmented into `nshoots` shooting segments stitched by Newton-iterated continuity, solved by `NonlinearSolveFirstOrder.NewtonRaphson` with `ForwardDiff` Jacobian and dense `LinearSolve.LUFactorization()`.

`MultipleShooting` dispatches on [`FreidlinWentzellHamiltonian`](@ref), not `CoupledSDEs`, because the shooting method operates on the LDT-derived deterministic Hamiltonian system rather than the original stochastic system. Passing a `CoupledSDEs` raises an `ArgumentError` pointing at the explicit construction:

```julia
H = FreidlinWentzellHamiltonian(sys)
res = minimize_geometric_action(H, x_initial, MultipleShooting(; nshoots = 10))
```

**Endpoint constraints.** Both endpoints must be *hyperbolic* fixed points of the drift (`‖b(x*)‖ ≈ 0`, `∂b(x*)` has no purely-imaginary eigenvalues). Non-fixed-point endpoints and non-hyperbolic fixed points are rejected with an `ArgumentError`. Use [`GeometricGradient`](@ref) (gMAM/sgMAM) in either regime.

**No interior fixed-point crossings.** The arclength reparametrization develops a singularity at any `(φ, p)` where `H_p = 0`, which on the `H = 0` shell is exactly `(x*, 0)` for `x*` a drift fixed point. The BVP-integrated portion of the path therefore must not pass through such a fixed point. The classical 1D bistable transition `xa = −1 → xb = +1` (which crosses the saddle `x = 0`) does not work as a single BVP; the user must split it into two `attractor → saddle` legs and sum the actions:

```julia
res_left  = minimize_geometric_action(H, x_init_left,  MultipleShooting())  # -1 → 0
res_right = minimize_geometric_action(H, x_init_right, MultipleShooting())  # +1 → 0
S_total = res_left.action + res_right.action
```

There is no runtime guard (cheap interior-fixed-point detection requires integrating, which the residual function does each Newton iteration anyway). If a user passes a through-saddle path, the BVP typically converges to a degenerate point with near-zero action; the symptom is `res.action` close to 0 when a nontrivial transition was expected.

**Warm-starting.** The BVP Newton iteration has a finite basin of attraction. Warm-starting from a `GeometricGradient` (gMAM/sgMAM) solve helps on nontrivial problems; `x_initial` follows the same `D × N` matrix convention as gMAM, so the swap is a one-line change.

The output `MinimumActionPath` carries the BVP-integrated path in `path`, the conjugate momentum in `generalized_momentum`, and the converged path length in `λ`.

```@docs
MultipleShooting
```

### `MinimumActionPath`
[(gMAM)](@ref "Geometric minimum action method (gMAM)") and [(sgMAM)](@ref "Simple geometric minimum action method (sgMAM)") return their output as a `MinimumActionPath` type:

```@docs
CriticalTransitions.MinimumActionPath
```

## Action functionals

### Freidlin-Wentzell action
```@docs
fw_action
```

### Geometric Freidlin-Wentzell action
```@docs
geometric_action
```

### Onsager-Machlup action
```@docs
om_action
```

For convenience, a general [`action`](@ref) function is available where the type of functional is set as an argument:

```@docs
action
```

## String method
The string method is a technique for finding transition paths between two states in a dynamical system [e_string_2007](@citet). The method represents the path as a "string" of points that connects the states and evolves it to minimize the drift along the path. The resulting tangent path is parallel to the drift of the system, i.e., the string method computes the heteroclinic orbit. For non-gradient systems (detailed -balance is broken), the heteroclinic orbit differs from the transition path, it does correctly predict, it correctly captures the deterministic dynamics from the saddle point onward ("downhill" portion of the path).
```@docs
string_method
```
