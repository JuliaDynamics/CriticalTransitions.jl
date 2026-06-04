# Large deviation theory

In the weak-noise limit, the probability of a stochastic system following a given path obeys a large deviation principle (LDP): transitions between metastable states concentrate, with overwhelming probability, on a single most probable path called the **instanton** or **maximum likelihood path**. The instanton minimizes a **rate functional** (the action), and its minimal value is the **quasipotential**, which sets the exponential transition rate and the most likely exit location.

This page follows that chain:

1. [The action](@ref "The action"): the rate functional being minimized.
2. [Finding the instanton](@ref "Finding the instanton"): minimizing the action to obtain the most probable path.
3. [The quasipotential](@ref "The quasipotential"): the minimized action as a field, and the transition rate.

For the underlying theory see [freidlin_random_1998](@citet) and [borner_climate_2025](@cite); for a survey of numerical methods see [grafke_numerical_2019](@citet).

!!! info "System class"
    The methods on this page apply to ``D``-dimensional stochastic dynamical systems of the form
    ```math
    \text{d} \mathbf{x} = \mathbf{b} (\mathbf{x}) \text{d}t + \sigma \mathbf{\Sigma} \text{d}\mathbf{W}_t \,,
    ```
    where the drift field ``\mathbf{b}`` may be non-gradient but the noise term must consist of Gaussian noise (``\mathbf{W}_t`` is a ``D``-dimensional vector of independent standard Wiener processes) and a constant covariance matrix ``\mathbf{Q} = \mathbf{\Sigma}\mathbf{\Sigma}^\top``. This is a special case of the broader class of noise types supported by [`CoupledSDEs`](@ref). The convention relating ``\sigma``, ``\mathbf{Q}`` and [`noise_strength`](@ref) is described under [Noise strength and the trace convention](@ref).

## The action

The action (or rate functional) assigns a cost to every candidate path; the instanton is its minimizer and the minimal cost is the leading-order ``-\log`` probability of the transition. Three forms are available:

- **Freidlin-Wentzell action**, [`fw_action`](@ref): the standard rate functional for a path traversed in a *fixed* transition time ``T``.
- **Geometric Freidlin-Wentzell action**, [`geometric_action`](@ref): a time-reparameterization-invariant form, used when ``T`` is not fixed.
- **Onsager-Machlup action**, [`om_action`](@ref): the finite-noise functional that adds the ``\mathcal{O}(\sigma^2)`` drift-divergence correction; it reduces to the Freidlin-Wentzell action as ``\sigma \to 0``. Additive noise only.

All three evaluate the integrand in the *canonical* covariance metric fixed by the [trace convention](@ref "Noise strength and the trace convention"). A practical consequence: [`fw_action`](@ref) and [`geometric_action`](@ref) are independent of the `noise_strength` chosen at construction, whereas [`om_action`](@ref) takes ``\sigma`` as an explicit argument because its correction scales with ``\sigma^2``; pass the ``\sigma`` at which you want ``-\log P[\phi]`` evaluated.

```@docs; canonical=false
fw_action
geometric_action
om_action
action
```

## Finding the instanton

The instanton is the minimizer of the action. Several methods have been proposed to compute it; in the weak-noise limit, the minimum action path coincides with the most probable transition path. While the minimum action method (MAM) is the most basic version, it is often beneficial to minimize the [geometric action](@ref "The action") via a time-independent version called gMAM. The problem can also be cast in a Hamiltonian form, implemented as simple gMAM (sgMAM), which can have numerical advantages. In gradient systems, minimum action paths between attractors coincide with heteroclinic orbits, which can be computed via the string method.

All three action-minimizers (MAM, gMAM, sgMAM) require only **autonomous** noise and reject **rank-deficient** diffusion; they all support additive, diagonal-multiplicative, and general-matrix multiplicative diffusion. The Onsager-Machlup functional (selectable in MAM via `functional = "OM"`) is the one exception: it is implemented only for additive noise and throws otherwise.

| Method | Use when |
|---|---|
| **MAM** ([`minimize_action`](@ref)) | The transition time ``T`` is fixed and you want the action at that ``T`` (or you specifically need the Onsager-Machlup functional, which only has a time-parameterized form). |
| **gMAM / sgMAM** ([`minimize_geometric_action`](@ref)) | ``T`` is unknown and you want a time-reparameterization-invariant instanton. Prefer **sgMAM** (Hamiltonian picture) when an analytic [`FreidlinWentzellHamiltonian`](@ref)`(H_x, H_p)` is available or when AD-based `jacobian(ds)` evaluations are cheap; prefer **gMAM** for the closer-to-textbook geometric-action formulation. |
| **Multiple shooting** ([`MultipleShooting`](@ref)) | Endpoints are hyperbolic fixed points and you want a high-accuracy boundary-value solution; best warm-started from a gMAM/sgMAM solve. |
| **String method** ([`string_method`](@ref)) | You want the deterministic heteroclinic orbit (typical use: **gradient** systems, where it coincides with the instanton). In non-gradient systems it generally differs from the most probable transition path. |

!!! info "Action minimization as an optimal control problem"
    All of the action minimizers below (MAM, gMAM, sgMAM) can equivalently be cast as **optimal control problems**: find a path ``\mathbf{x}`` and a control ``\mathbf{u}`` that minimize a Freidlin-Wentzell-type action
    ```math
    \int \| \mathbf{u}(t) - \mathbf{b}(\mathbf{x}(t)) \|^2 \, \text{d}t
    ```
    subject to the path dynamics ``\dot{\mathbf{x}}(t) = \mathbf{u}(t)`` and the boundary conditions ``\mathbf{x}(0) = \mathbf{x}_A``, ``\mathbf{x}(T) = \mathbf{x}_B`` (with ``T`` either fixed, as in MAM, or eliminated by reparameterization, as in gMAM/sgMAM). In this form the problem can be solved with dedicated optimal-control packages such as [`OptimalControl.jl`](https://control-toolbox.org/OptimalControl.jl); see the [control-toolbox.org MAM tutorial](https://control-toolbox.org/Tutorials.jl/stable/tutorial-mam.html) for a worked example on the Maier-Stein system.

#### Variants and extensions
The literature contains a number of extensions of MAM-type methods that may be relevant depending on the model class and numerical difficulties. These variants are not currently implemented in `CriticalTransitions.jl`, but serve as useful pointers:

- **tMAM / optimal linear time scaling**: avoids explicit optimization over the transition time by introducing an optimal linear time scaling; can be combined with adaptivity in time discretization [wan_tmam_2015](@citet).
- **Adaptive MAM**: uses a moving-mesh strategy to concentrate grid points in dynamically important portions of the path, improving efficiency and robustness [zhou_adaptive_mam_2008](@citet).
- **Non-Gaussian (jump / Lévy) noise**: for systems driven by jump noise, the rate function and path optimization problem differ from the Freidlin-Wentzell diffusive setting; see e.g. an optimal-control-based approach in [wei_most_likely_jumps_2023](@citet).
- **Multiplicative / state-dependent noise**: state-dependent diagonal and full-matrix multiplicative noise is supported by both gMAM and sgMAM. The diffusion tensor ``a(x)`` is classified once when the sgMAM cache is built and the resulting inner loop dispatches at compile time on the cache type (see [Developer notes & internals](@ref)); the Hamiltonian formulation is given in [grafke_small_random_2017](@citet). Rank-deficient (degenerate) noise is rejected at cache build; see [grafke_small_random_2017](@citet) for the rank-deficient extension that would be needed to support it.

### Minimum action method (MAM)
Minimization of the specified action functional using the optimization algorithm of `Optimization.jl`. See also [e_minimum_2004](@citet).

```@docs; canonical=false
minimize_action
```

### Geometric minimum action method (gMAM / sgMAM)
gMAM minimizes the geometric action following [heymann_pathways_2008, heymann_geometric_2008](@citet). It reformulates MAM to avoid the double optimization over both the action and the transition time, using a [geometric action](@ref "The action") functional that is independent of the time parametrization of the path. This reparameterization invariance makes the method more robust and computationally efficient, particularly for multiscale systems.

The simple gMAM (sgMAM) following [grafke_long_2017](@citet) reduces the complexity of the original gMAM by requiring only first-order derivatives of the underlying Hamiltonian optimization formulation. This simplifies the numerical treatment and computational complexity. The implementation performs a constrained gradient descent on the Hamiltonian system, supporting autonomous diffusions with additive, diagonal-multiplicative, or general-matrix multiplicative noise.

The gradient step size is set by the [`GeometricGradient`](@ref) optimizer passed to [`minimize_geometric_action`](@ref); for difficult problems the probe-based [`AdaptiveGeometricGradient`](@ref) variant is more robust (see the adaptive step-size example).

#### Hamiltonian picture

The Freidlin-Wentzell rate function admits an equivalent Hamiltonian formulation with conjugate momentum ``p``:

```math
H(\varphi, p) = \langle b(\varphi), p \rangle + \tfrac{1}{2}\, \langle p,\, a(\varphi)\, p \rangle,
```

where ``a(x) = \sigma(x)\sigma(x)^{\top}``. The action along an instanton (where ``H \equiv 0``) reduces to ``S[\varphi] = \int_0^T \langle p, \mathrm{d}\varphi\rangle``, the form used by [`minimize_geometric_action`](@ref) when the input is a [`FreidlinWentzellHamiltonian`](@ref). Rank-deficient ``a(x)`` is rejected when the per-path cache is built. The compile-time classification of ``a(x)`` is described under [Developer notes & internals](@ref).

```@docs; canonical=false
minimize_geometric_action
FreidlinWentzellHamiltonian
```

#### Performance notes (sgMAM)

sgMAM repeatedly evaluates `H_p(x, p)` and `H_x(x, p)` along a discretized path. If this allocates, sgMAM will be slow. The key reason for performance differences:

- `FreidlinWentzellHamiltonian(ds::CoupledSDEs)` typically relies on `jacobian(ds)` (often automatic differentiation unless you provide an analytic Jacobian) and evaluates it pointwise along the path.
- A hardcoded `FreidlinWentzellHamiltonian(H_x, H_p)` with analytic expressions operating on the full `D×Nt` path matrix usually allocates far less.

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

The same "vectorized + allocation-free inner loop" principle also tends to make [`string_method`](@ref) faster when used with `FreidlinWentzellHamiltonian`.

### Multiple shooting

The [`MultipleShooting`](@ref) optimizer treats the instanton as a boundary value problem on arclength-reparameterized Hamilton equations
```math
\frac{\mathrm{d}\varphi}{\mathrm{d}s} = \alpha\,H_p(\varphi, p),\qquad
\frac{\mathrm{d}p}{\mathrm{d}s} = -\alpha\,H_x(\varphi, p),\qquad
\alpha = L\,/\,\|H_p\|
```
on `s ∈ [0, 1]`, with the total path length `L` a Newton unknown so that `‖dφ/ds‖ ≡ L` pointwise. The BVP is segmented into `nshoots` shooting segments stitched by Newton-iterated continuity; the solver internals are described under [Developer notes & internals](@ref).

`MultipleShooting` dispatches on [`FreidlinWentzellHamiltonian`](@ref), not `CoupledSDEs`, because the shooting method operates on the LDT-derived deterministic Hamiltonian system rather than the original stochastic system. Passing a `CoupledSDEs` raises an `ArgumentError` pointing at the explicit construction:

```julia
H = FreidlinWentzellHamiltonian(sys)
res = minimize_geometric_action(H, x_initial, MultipleShooting(; nshoots = 10))
```

**Endpoint constraints.** Both endpoints must be *hyperbolic* fixed points of the drift (`‖b(x*)‖ ≈ 0`, `∂b(x*)` has no purely-imaginary eigenvalues). Non-fixed-point endpoints and non-hyperbolic fixed points are rejected with an `ArgumentError`. Use [`GeometricGradient`](@ref) (gMAM/sgMAM) in either regime.

**No interior fixed-point crossings.** The arclength reparameterization develops a singularity at any `(φ, p)` where `H_p = 0`, which on the `H = 0` shell is exactly `(x*, 0)` for `x*` a drift fixed point. The BVP-integrated portion of the path therefore must not pass through such a fixed point. The classical 1D bistable transition `xa = −1 → xb = +1` (which crosses the saddle `x = 0`) does not work as a single BVP; the user must split it into the noise-driven `attractor → saddle` leg and the deterministic `saddle → attractor` leg. Only the uphill `attractor → saddle` leg carries action; the downhill `saddle → attractor` leg follows the deterministic drift and contributes zero Freidlin-Wentzell action. The transition action is therefore the action of the single uphill leg (and, when several competing barrier legs exist, the *maximum* over them), not the sum of both `attractor → saddle` solves:

```julia
res_left  = minimize_geometric_action(H, x_init_left,  MultipleShooting())  # -1 → 0 (uphill)
res_right = minimize_geometric_action(H, x_init_right, MultipleShooting())  # +1 → 0 (uphill)
# -1 → +1 transition: uphill -1 → 0, then deterministic 0 → +1 (zero action).
S_transition = res_left.action
# Symmetric case: res_left.action ≈ res_right.action.
```

There is no runtime guard for this. If a user passes a through-saddle path, the BVP typically converges to a degenerate point with near-zero action; the symptom is `res.action` close to 0 when a nontrivial transition was expected.

**Warm-starting.** The BVP Newton iteration has a finite basin of attraction. Warm-starting from a `GeometricGradient` (gMAM/sgMAM) solve helps on nontrivial problems; `x_initial` follows the same `D × N` matrix convention as gMAM, so the swap is a one-line change.

The output [`MinimumActionPath`](@ref CriticalTransitions.MinimumActionPath) carries the BVP-integrated path in `path`, the conjugate momentum in `generalized_momentum`, and the converged path length in `λ`.

```@docs; canonical=false
MultipleShooting
```

### String method
The string method finds transition paths between two states by representing the path as a "string" of points connecting the states and evolving it to minimize the drift along the path [e_string_2007](@citet). The resulting tangent path is parallel to the drift of the system, i.e., the string method computes the heteroclinic orbit. For non-gradient systems (where detailed balance is broken) the heteroclinic orbit differs from the transition path, but it correctly captures the deterministic dynamics from the saddle point onward (the "downhill" portion of the path).

```@docs; canonical=false
string_method
```

### `MinimumActionPath`
[gMAM/sgMAM](@ref "Geometric minimum action method (gMAM / sgMAM)") and [multiple shooting](@ref "Multiple shooting") return their output as a `MinimumActionPath`:

```@docs; canonical=false
CriticalTransitions.MinimumActionPath
```

## The quasipotential

The quasipotential ``U_A(\mathbf{x})`` is the minimal action required to reach a state ``\mathbf{x}`` from an attractor ``\mathbf{x}_A``. It generalizes the potential of gradient systems to non-gradient drift and sets both the exponential transition rate (``\sim e^{-U_A/\sigma^2}``) and the most likely exit location. There are two ways to obtain it:

- **Along a single instanton**: minimize the action between two states with one of the [minimizers above](@ref "Finding the instanton") and read off the action value.
- **As a field**: the Ordered Line Integral Method (OLIM) [Dahiya2018OLIM](@cite) computes the entire quasipotential field on a Cartesian grid in a single Dijkstra-style sweep, given a stable attractor ``\mathbf{x}_A``. The dimension ``D`` is taken from the `CoupledSDEs` type parameter; the method is most accurate for ``D = 2, 3``.

```julia
using CriticalTransitions, StaticArrays
f(x, p, t) = SVector(x[1] - x[1]^3 - 5 * x[1] * x[2]^2,
                      -(1 + x[1]^2) * x[2])
sys  = CoupledSDEs(f, [-1.0, 0.0]; noise_strength = 0.3)
grid = CartesianGrid((-1.5, 1.5, 121), (-1.0, 1.0, 81))
qp   = quasipotential(sys, grid, [-1.0, 0.0])
qp.U[CartesianIndex(61, 41)]   # U at the saddle (0, 0)
```

See `examples/quasipotential_maierstein.jl` for a worked example with contour plot.

```@docs; canonical=false
quasipotential
CriticalTransitions.QuasiPotential
CriticalTransitions.BackRef
```
