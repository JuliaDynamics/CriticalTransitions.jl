# Multiple-shooting BVP for the geometric Freidlin-Wentzell instanton

**Status:** Draft, awaiting user review.
**Date:** 2026-05-22.
**Issues addressed:** [#195](https://github.com/JuliaDynamics/CriticalTransitions.jl/issues/195) (Action plot method, repurposed to mean "multiple-shooting BVP for the FW instanton"; the forward-shooting / action-plot variant of Beri 2005 is deferred to a follow-up).
**Depends on:** [`2026-05-22-fw-hamiltonian-and-multiplicative-noise-design.md`](2026-05-22-fw-hamiltonian-and-multiplicative-noise-design.md). The shooting method consumes `FreidlinWentzellHamiltonian{IIP, D, Hx, Hp, AF, NS}` from that spec and cannot merge before it.

## 1. Goals and non-goals

### Goals

* Add a `MultipleShooting <: AbstractActionOptimizer` optimizer for
  `minimize_geometric_action(::FreidlinWentzellHamiltonian, x_initial, MultipleShooting(...); kwargs...)`.
* Formulate the BVP on the **geometric (arclength) Hamilton equations**. The path parameter is arclength `s ∈ [0, 1]`. The terminal time `T` does not appear in the user-facing API and is not an internal Newton unknown.
* Auto-classify endpoints by `‖b(x)‖`. Fixed-point endpoints (attractors or saddles) use **linearized endpoint segments** built from the eigendecomposition of the linearized Hamiltonian at the fixed point. Free endpoints (non-zero drift) use direct Dirichlet boundary conditions.
* Support `AdditiveNoise`, `DiagonalNoise`, and `GeneralNoise` uniformly. The shooting code calls `H_x`, `H_p` from the FW Hamiltonian as black boxes, so multiplicative-noise handling is inherited.
* Compile-time dispatch on the noise shape: `MultipleShootingWorkspace{D, NS}` carries `NS` as a type parameter, mirroring the FW spec.
* Output is the existing `MinimumActionPath` struct, extended with BVP-specific diagnostic metadata (residual norm, nlsolve iteration count, `H = 0` invariant violation).

### Non-goals

* `MultipleShooting` on `::CoupledSDEs` directly. The user constructs `FreidlinWentzellHamiltonian(sys)` explicitly. A custom `Base.showerror` hint guides the migration. The conceptual point (shooting acts on the LDT Hamiltonian, not the original SDE) is encoded in the dispatch table.
* Time-domain (`T` as Newton unknown) shooting. Out of scope. Arclength is the only formulation we ship.
* Forward shooting from the unstable manifold (Beri et al. 2005 action plot). Separate feature, separate PR.
* Floquet-extended (periodically driven, non-autonomous) shooting (arXiv:2203.14329). Out of scope. The autonomous version is a prerequisite.
* AD for Hamiltonian Jacobians. Defer. Use `BoundaryValueDiffEqShooting`'s default FD-based Jacobian construction. A future PR can wire in `DifferentiationInterface.jl` via a kwarg without a breaking change.
* `string_method`, `minimize_action`, `om_action`, `fw_action`: untouched.

## 2. Architecture

### 2.1 `MultipleShooting` optimizer type

```julia
struct MultipleShooting{S, J} <: AbstractActionOptimizer
    nshoots::Int                # number of shooting segments (default 10)
    ode_solver::S               # default Tsit5()
    nlsolve::J                  # default nothing => BoundaryValueDiffEqShooting picks
    abstol::Float64             # BVP solver abstol (default 1e-8)
    reltol::Float64             # BVP solver reltol (default 1e-6)
end

MultipleShooting(; nshoots = 10, ode_solver = Tsit5(),
                   nlsolve = nothing, abstol = 1e-8, reltol = 1e-6) = ...
```

`nshoots` is the BVP segment count, not the output path resolution. The output path resolution `N` is derived from `size(x_initial, 2)` (the user-supplied initial guess matrix), exactly as the existing `minimize_geometric_action(::CoupledSDEs, x_initial, ::GeometricGradient)` does.

### 2.2 Endpoint classification

Computed once for each of `xa = x_initial[:, 1]`, `xb = x_initial[:, end]` at the start of the call:

```julia
function _classify_endpoint(H::FreidlinWentzellHamiltonian{IIP, D},
                            x::AbstractVector) where {IIP, D}
    b = H.H_p(x, zeros(D))                       # drift, p = 0
    norm(b) < sqrt(eps(eltype(x))) ? :fixed_point : :free
end
```

The tolerance follows the noise-shape classifier of the FW spec for consistency. Both endpoints are classified independently. The four combinations (`{fp, free} × {fp, free}`) are all supported; the linearized-segment branch is taken per-side.

### 2.3 Output

The BVP solution (a `BVPSolution` from `BoundaryValueDiffEqShooting`) is sampled at `N` uniform arclength nodes via its built-in interpolant, then packaged into the existing `MinimumActionPath`:

* `path :: Matrix{T}`, size `D × N`. Configuration-space samples only; `p` is internal.
* `action :: T`, computed via `geometric_action(b, path, arclength; A = inv∘a)` on the sampled path (Section 3.4).
* `arclength :: Vector{T}`, uniform on `[0, 1]`.
* `metadata :: NamedTuple` (existing field extended), containing:
  * `bvp_residual_norm :: T`,
  * `nlsolve_iterations :: Int`,
  * `H_invariant_max :: T` (`max_s |H(x(s), p(s))|`).

Existing consumers of `MinimumActionPath` are unaffected; metadata gains optional fields keyed by symbol.

### 2.4 Geometric arclength Hamilton equations

The ODE integrated between shooting nodes:

```
dx/ds = λ(s) · H_p(x, p)
dp/ds = −λ(s) · H_x(x, p)
λ(s) = ‖H_p(x, 0)‖_{a⁻¹} / ‖H_p(x, p)‖_{a⁻¹}
```

where `‖·‖_{a⁻¹}` denotes the `a(x)⁻¹`-weighted norm. The `λ` closure is derived from the on-shell condition `H(x, p) = 0` and removes the need for a DAE formulation. `H_x`, `H_p`, `a` come from `FreidlinWentzellHamiltonian` and already encode the noise shape (additive, diagonal, general). The shooting code does not branch on `NS` inside the ODE RHS.

For `AdditiveNoise` the closure simplifies (`‖H_p(x, 0)‖_{a⁻¹} = ‖b‖_{a⁻¹}` is constant in `s`) and `λ(s)` reduces to a single division.

## 3. Algorithm internals

### 3.1 Linearized endpoint segments

For each fixed-point endpoint `x*`, replace the corresponding shooting segment with the analytical solution of the linearized Hamiltonian flow around `(x*, 0)`. This avoids the `0/0` degeneracy in `λ(s)` at the endpoint.

Define

```
J  = ∇b(x*)            (drift Jacobian; FD around x* with step h = √eps)
A* = a(x*)             (noise covariance at the fixed point)
M  = [ J         A*  ]
     [ 0       −Jᵀ   ]
```

`M` is the Jacobian of the time-domain Hamilton equations at `(x*, 0)`. It has Hamiltonian structure: eigenvalues come in `±μ` pairs. Assuming `x*` is hyperbolic for the drift (`J` has no purely imaginary eigenvalues), `M` has exactly `D` eigenvalues with `Re(μ) > 0` (unstable) and `D` with `Re(μ) < 0` (stable).

Compute eigendecomposition once per fixed-point endpoint at call time:

```
M · U_u = U_u · Λ_u     # D unstable eigenpairs
M · U_s = U_s · Λ_s     # D stable eigenpairs
```

**Outgoing instanton at `xa = x*_a`:**

```
[x(s) − xa; p(s)] = U_u_a · exp(Λ_u_a · (τ(s) − τ_max)) · c_a,    s ∈ [0, s_1]
```

with `c_a ∈ ℝ^D` the unstable-mode coefficients (Newton unknowns), `U_u_a` the matrix of unstable eigenvectors of `M_a`, and `τ(s)` a monotone arclength-to-linearization-time map. `τ(s_1) = τ_max` (segment hand-off to the nonlinear BVP integration) and `τ(0) → −∞` in the formal limit (the fixed point sits at infinite linearized "time"). In practice the implementation picks a finite linearization horizon `τ_min` so that `‖[x(0); p(0)] − [xa; 0]‖ < ε_lin` (default `ε_lin = 1e-8`). Both `τ_min` and the precise form of `τ(s)` are implementation choices; a sensible default is "uniform arclength sampling on the linearized curve", which is a 1D root-find solved once per Newton step.

At `s = 0`, the parameterization gives `x ≈ xa` to within `ε_lin`. The residual contribution at `s = 0` is therefore implicit (absorbed into the parameterization) rather than an explicit Newton equation.

**Incoming instanton at `xb = x*_b`:**

```
[x(s) − xb; p(s)] = U_s_b · exp(Λ_s_b · (τ(s) − τ_min_b)) · c_b,    s ∈ [s_{N_seg − 1}, 1]
```

symmetric. The stable eigenvalues have `Re < 0`, so the exponential decays into `(xb, 0)` as `τ(s) → +∞`. Same finite-horizon trick: `τ(1)` is set to a large positive value so the path arrives within `ε_lin` of `(xb, 0)`.

Free endpoints get no linearized segment; `(x_0, p_0)` (or `(x_N, p_N)`) are direct Newton unknowns with the boundary residual `x_0 = xa` (or `x_N = xb`) replacing the matching residual.

Non-hyperbolic fixed point (any `Re(μ) ≈ 0` of `J`) raises an `ArgumentError` at construction:

```
"linearized endpoint segment requires hyperbolic fixed point at xa (or xb); 
 either perturb the endpoint slightly off the bifurcation or use 
 GeometricGradient (gMAM/sgMAM) which handles non-hyperbolic endpoints via discretization"
```

### 3.2 BVP unknowns and residuals

With `N_seg = nshoots`, arclength grid `s_0 = 0 < s_1 < ... < s_{N_seg} = 1` uniform. The node states are `y_i = (x_i, p_i) ∈ ℝ^{2D}` for `i = 0, 1, ..., N_seg`. There are `N_seg + 1` nodes and `N_seg` segments.

**Unknowns by node:**

| Node            | Endpoint kind     | Unknown count | Parameterization                            |
|-----------------|-------------------|---------------|---------------------------------------------|
| `y_0`           | Fixed point       | `D` (`c_a`)   | linearized: `y_0 = U_u_a · exp(Λ_u_a · (τ_min − τ_max)) · c_a + [xa; 0]` |
| `y_0`           | Free              | `D` (`p_0`)   | direct: `y_0 = (xa, p_0)`                   |
| `y_i`, `i = 1..N_seg−1` | (interior) | `2D`          | direct: `y_i = (x_i, p_i)`                  |
| `y_{N_seg}`     | Fixed point       | `D` (`c_b`)   | linearized: symmetric to `y_0`              |
| `y_{N_seg}`     | Free              | `D` (`p_{N_seg}`) | direct: `y_{N_seg} = (xb, p_{N_seg})`   |

Total unknowns: `D + 2D · (N_seg − 1) + D = 2D · N_seg`, independent of endpoint kinds.

**Residuals:** all boundary conditions are absorbed into the node parameterization (free endpoints fix `x_0 = xa`, `x_{N_seg} = xb` by setting them as constants; fixed-point endpoints fix `x ≈ x*` to within `ε_lin` by the linearized parameterization). No explicit Dirichlet boundary residuals appear.

The remaining residuals are `N_seg` segment-matching equations, one per segment:

```
ϕ_i(y_i, y_{i+1}) := integrate_geometric_Hamilton_ODE(y_i, s_i → s_{i+1}) − y_{i+1} = 0,
                                                                  i = 0, 1, ..., N_seg − 1
```

Each `ϕ_i` returns `2D` residuals. Total: `2D · N_seg`.

Residual count equals unknown count by construction. The Newton system is square. The BVP is well-posed when the instanton between `xa` and `xb` exists and is locally unique (the generic case for hyperbolic fixed points and well-chosen free endpoints; degenerate cases such as bifurcations of the fixed point are out of scope).

### 3.3 Per-shape compile-time dispatch

```julia
struct MultipleShootingWorkspace{IIP, D, NS<:NoiseShape, HT, BVP, CB, CC}
    H::HT                       # FreidlinWentzellHamiltonian{IIP, D, ..., NS}
    nshoots::Int
    endpoint_a::Symbol          # :fixed_point or :free
    endpoint_b::Symbol
    linearization_a::CB         # nothing or (U_u, Λ_u, c_a_init); for fixed-point xa
    linearization_b::CB         # nothing or (U_s, Λ_s, c_b_init); for fixed-point xb
    bvp_problem::BVP            # BVProblem from BoundaryValueDiffEqShooting
    bvp_cache::CC               # LinearSolve / NonlinearSolve cache reused across iters
end
```

`NS` flows in from `H` and parameterizes the workspace. The ODE RHS closure captures `H.H_x` and `H.H_p`, which themselves dispatch on `NS` per the FW spec. No runtime `if NS isa ...` branches inside the integrator hot loop. Verified by `@code_warntype` on the inner RHS evaluation.

### 3.4 Action computation

After BVP convergence:

```julia
path = sample_bvp_at_arclength_nodes(sol, N)   # D × N matrix
S = geometric_action(b, path, arclength; A = inv∘a)
```

This reuses the existing `geometric_action(b, path, arclength; A)` integrand from the FW spec section 4.2. The trace-normalized covariance convention and PR #329 action canonicalization are inherited automatically; no parallel implementation to keep in sync.

We deliberately do **not** accumulate `S` along the BVP ODE integration. The BVP loop is for path convergence only; action is a measurement on the converged path.

### 3.5 `H = 0` invariant check

After BVP convergence, evaluate

```
H_max = maximum_s |H(x(s), p(s))|
```

at, say, `4 N_seg` sample points along the BVP interpolant. If `H_max > 1e-6` (configurable via `MultipleShooting.invariant_tol`, default `1e-6`), emit a `@warn` with the maximum violation and the call-site information. The warning does not abort: a small violation typically reflects BVP tolerance, not algorithmic error.

For the linearized endpoint segments, `H(x(s), p(s)) ≡ 0` is exact on the linearization (by construction; the linearized Hamiltonian shares its zero level set with the unstable/stable invariant manifolds). So `H_max` only probes the BVP-integrated interior.

### 3.6 Performance commitments

In scope for this PR, mirroring FW spec section 3.4:

* **Type stability**: `minimize_geometric_action(::FreidlinWentzellHamiltonian, ..., ::MultipleShooting)` is fully type-stable end-to-end. Verified via `Test.@inferred` on one path per `NS` branch in the test suite.
* **`MultipleShootingWorkspace` is fully concrete**: all type parameters resolved at construction.
* **`bvp_cache` reused across BVP iterations**: `BoundaryValueDiffEqShooting` already handles this internally for the Jacobian/linear solve; we configure `nlsolve = nothing` to let the default cache reuse kick in. No `@allocations` regression target for the inner Newton loop because that lives in `BoundaryValueDiffEqShooting` and is its responsibility.
* **No multi-allocation regression in the ODE RHS**: `λ(s)` closure and `H_p`, `H_x` calls are all allocating-free for `AdditiveNoise` and `DiagonalNoise` (verified by `@allocations` test on the RHS at a fixed `(x, p)`).
* **Benchmark**: new `benchmarks/shooting_maierstein.jl` covering additive Maier-Stein. Threshold: shooting wall time within `2×` of gMAM on the same problem at `nshoots = 10`, `N = 100`. The factor reflects that shooting does a Newton solve while gMAM does fixed-step descent; we are not claiming shooting is faster, only that it's competitive enough to be useful.

## 4. API surface

### 4.1 Dispatch table

| `sys` type | `optimizer` type | Algorithm | Status |
|------------|------------------|-----------|--------|
| `CoupledSDEs` | `GeometricGradient` / `Adam` | gMAM | existing |
| `FreidlinWentzellHamiltonian` | `GeometricGradient` | sgMAM | from FW spec |
| `FreidlinWentzellHamiltonian` | **`MultipleShooting`** | **shooting (this PR)** | new |
| `CoupledSDEs` | `MultipleShooting` | (rejected, `MethodError`) | new |

### 4.2 `MethodError` hint for `CoupledSDEs`

Custom error path so users hit `MultipleShooting` on a `CoupledSDEs` see:

```
ERROR: MultipleShooting acts on FreidlinWentzellHamiltonian, not CoupledSDEs.
The shooting method operates on the LDT-derived deterministic Hamiltonian system,
not the original stochastic system. Convert explicitly:

    H = FreidlinWentzellHamiltonian(sys)
    minimize_geometric_action(H, x_initial, MultipleShooting(...))
```

Implementation: a small `Base.showerror` extension following the pattern the FW spec uses for closing #263.

### 4.3 Exports

`src/CriticalTransitions.jl` adds:

* `MultipleShooting`.

All other exports unchanged. No removals from this PR.

## 5. Testing strategy

### 5.1 File: `test/largedeviations/multiple_shooting.jl`

**Cross-validation (Section 4.1 of the brainstorm):**

* Maier-Stein `β = 1`, attractor `→` saddle `→` attractor, additive isotropic. Shooting action equals analytic `S = 1/2` to `rtol = 1e-3`. Path matches gMAM to `rtol = 1e-2` after re-interpolation onto a common grid.
* Maier-Stein attractor `→` saddle (half instanton): action equals `1/4`. Endpoint classifier: both `:fixed_point`.
* Maier-Stein free `→` free: `xa`, `xb` chosen off any fixed point. Shooting vs gMAM, `rtol = 1e-2`. Endpoint classifier: both `:free`.
* 1D OU additive: analytic instanton; action and path match to `rtol = 1e-4`.

**Multiplicative noise:**

* 1D OU multiplicative `DiagonalNoise`: shooting vs analytic Simpson, `rtol = 1e-3`. Shooting vs gMAM, `rtol = 0.05`.
* 2D off-diagonal `σ(x)` `GeneralNoise`: shooting vs gMAM, `rtol = 0.1`.
* Rotated Maier-Stein (constant non-diagonal `Σ`, classified `GeneralNoise`): action matches diagonalized form to `rtol = 1e-3`.

**Invariants:**

* `H = 0` along solution: `max_s |H(x(s), p(s))| < 1e-6` on each converged case. Asserted once per converged test.
* Convergence in `nshoots`: action stable across `nshoots ∈ {5, 10, 20}` to `rtol = 1e-3`.
* Reverse-direction symmetry: action `xa → xb` equals action `xb → xa` to `rtol = 1e-4` (when both endpoints are fixed points).

**Compile-time dispatch:**

* `MultipleShootingWorkspace{D, NS}` carries the right `NS` for each test SDE; introspect via `typeof`.
* `Test.@inferred` on the entry-point call for one path per `NS` branch.

### 5.2 File: `test/largedeviations/multiple_shooting_rejection.jl`

* `minimize_geometric_action(::CoupledSDEs, x, MultipleShooting(...))` raises `MethodError`. The error message contains the strings `"FreidlinWentzellHamiltonian"` and `"shooting"`.
* Non-hyperbolic fixed point (constructed system with a zero eigenvalue of `∇b`): `MultipleShooting` raises `ArgumentError` mentioning `"hyperbolic"`.

### 5.3 Coverage matrix

| Endpoint pair | Noise | Cross-checked against |
|---|---|---|
| attractor `↔` attractor (via saddle) | additive iso | analytic, gMAM, sgMAM |
| attractor `→` saddle | additive iso | analytic (half action) |
| free `↔` free | additive iso | gMAM |
| attractor `↔` attractor | diagonal multiplicative | analytic Simpson, gMAM |
| attractor `↔` attractor | general multiplicative | gMAM |
| attractor `↔` attractor | rotated `Σ` (constant non-diagonal) | diagonalized-basis form |

## 6. File-level change list

### New files

* `src/largedeviations/multiple_shooting.jl`: `MultipleShooting` struct, `MultipleShootingWorkspace`, `_classify_endpoint`, `_build_linearization`, geometric-Hamilton ODE RHS, BVP assembly via `BoundaryValueDiffEqShooting.MultipleShooting`, action computation wrapper, `H`-invariant check.
* `test/largedeviations/multiple_shooting.jl`: cross-validation, multiplicative, invariants, compile-time dispatch (Section 5.1).
* `test/largedeviations/multiple_shooting_rejection.jl`: `MethodError` hint and non-hyperbolic rejection (Section 5.2).
* `examples/shooting_Maierstein.jl`: worked example. Builds `FreidlinWentzellHamiltonian(sys)` for the Maier-Stein SDE, runs both `GeometricGradient` and `MultipleShooting`, plots both paths and prints both actions.
* `benchmarks/shooting_maierstein.jl`: timing benchmark on additive Maier-Stein. Hooks into the existing `benchmark/` driver.

### Modified files

* `src/largedeviations/minimize_geometric_action.jl`: add the `(::FreidlinWentzellHamiltonian, x_initial, ::MultipleShooting)` dispatch (delegates to the constructor in `multiple_shooting.jl`). Add the `(::CoupledSDEs, x_initial, ::MultipleShooting)` `MethodError` hint via `Base.showerror`.
* `src/largedeviations/methods.jl`: declare `MultipleShooting <: AbstractActionOptimizer` (or equivalent supertype). Confirm at implementation time which is the existing common base.
* `src/largedeviations/MinimumActionPath.jl`: extend the metadata `NamedTuple` to optionally carry `bvp_residual_norm`, `nlsolve_iterations`, `H_invariant_max`. Strictly additive; existing consumers unaffected.
* `src/CriticalTransitions.jl`: add `MultipleShooting` to the export list.
* `Project.toml`: add `BoundaryValueDiffEqShooting` dependency. We depend on the shooting-only sub-package rather than the umbrella `BoundaryValueDiffEq` because we only use `MultipleShooting` (and possibly `Shooting`). Skipping the MIRK / FIRK / Ascher sub-packages keeps precompilation lighter. Compat bound matches the version currently used alongside `OrdinaryDiffEq` in the SciML stack (verified at implementation time). Same ecosystem as the existing `LinearSolve` / `OrdinaryDiffEq` deps; no new external system pulled in.
* `test/runtests.jl`: include `multiple_shooting.jl` and `multiple_shooting_rejection.jl`.
* `docs/src/man/largedeviations.md`: new "Multiple shooting" subsection. Walks through the FW Hamiltonian construction, the `MultipleShooting(...)` optimizer, endpoint auto-classification, and a worked Maier-Stein example. Explicitly documents *why* the user must construct `FreidlinWentzellHamiltonian` (the conceptual point baked into the dispatch).

### Untouched

* `src/largedeviations/sgMAM.jl`, `src/largedeviations/action.jl`, `src/largedeviations/string_method.jl`, `src/largedeviations/utils.jl`: no changes. This PR is purely additive at the algorithm layer.
* gMAM behavior on `CoupledSDEs` and on `FreidlinWentzellHamiltonian`: untouched.

## 7. Risks and migration

### Breaking changes

None for existing API. `MultipleShooting` is purely additive.

### Convergence risks

* **Poor initial guess**: BVP Newton convergence depends on the initial path being in the basin of attraction of the instanton. Mitigation: documentation recommends warm-starting from a gMAM solution for hard problems. The `x_initial` API matches gMAM's so this is a one-line swap for the user.
* **Multiple instantons**: if the system has multiple instantons between `xa` and `xb`, BVP shooting converges to whichever is nearest the initial guess. Same caveat as gMAM/sgMAM. Documentation flags this.
* **Non-hyperbolic endpoints**: explicit `ArgumentError` rejection (Section 3.1). Users are pointed at `GeometricGradient` as the fallback.

### Performance risks

* **Jacobian build cost**: `BoundaryValueDiffEqShooting` defaults to FD Jacobian, which for large `D · N_seg` can dominate. The benchmark threshold (`2×` gMAM) is generous to account for this. If the FD Jacobian is the bottleneck, the follow-up is AD via `DifferentiationInterface.jl` (configurable kwarg, not a breaking change).
* **Type instability**: would be a release blocker. Covered by `Test.@inferred` in the test suite.

### Out-of-scope follow-ups

* **Forward shooting / action plot (Beri 2005)**: separate PR. Reuses the geometric ODE RHS and the FW Hamiltonian; replaces the BVP solver with a sweep over initial conditions on the unstable manifold.
* **Floquet-extended shooting** for periodically driven systems (arXiv:2203.14329): builds on action plot, needs the autonomous version first.
* **AD-based Jacobian** via `DifferentiationInterface.jl`: drop-in once a workload demands it.
* **Continuation in a parameter**: warm-start each parameter step from the previous BVP solution. Cheap follow-up if users ask for it.

## 8. References

* Heymann, M. and Vanden-Eijnden, E. (2008). *The geometric minimum action method: a least action principle on the space of curves*. Communications on Pure and Applied Mathematics 61(8). Geometric action formulation, eqs 2.4 and 3.6.
* Grafke, T., Schäfer, T., Vanden-Eijnden, E. (2017). *Long term effects of small random perturbations*. gMAM as a gradient descent on the geometric action.
* Grafke, T. and Vanden-Eijnden, E. (2019). *Numerical computation of rare events via large deviation theory*. Chaos 29 063118. [arXiv:1812.00681](https://arxiv.org/abs/1812.00681). Hamiltonian formulation of the FW instanton, eqs 17 and 28.
* Beri, S., Mannella, R., Luchinsky, D. G., Silchenko, A. N., McClintock, P. V. E. (2005). *Solution of the boundary value problem for optimal escape in continuous stochastic systems and maps*. Forward shooting / action plot method.
* Lindner, M., Hellman, F., Kurths, J. (2018). *Asymptotic instanton calculus for stochastic differential equations*. Multiple-shooting solvers for LDT instantons.
* Internal: FW Hamiltonian + multiplicative-noise spec [`2026-05-22-fw-hamiltonian-and-multiplicative-noise-design.md`](2026-05-22-fw-hamiltonian-and-multiplicative-noise-design.md), PR #329 (action conventions), issue #195.
