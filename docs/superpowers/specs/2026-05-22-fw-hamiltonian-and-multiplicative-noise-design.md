# Freidlin-Wentzell Hamiltonian, multiplicative noise, unified entry point

**Status:** Draft, awaiting user review.
**Date:** 2026-05-22.
**Issues addressed:** [#263](https://github.com/JuliaDynamics/CriticalTransitions.jl/issues/263), [#275](https://github.com/JuliaDynamics/CriticalTransitions.jl/issues/275), [#326](https://github.com/JuliaDynamics/CriticalTransitions.jl/issues/326), [#328](https://github.com/JuliaDynamics/CriticalTransitions.jl/issues/328).
**Punted to follow-up:** [#325](https://github.com/JuliaDynamics/CriticalTransitions.jl/issues/325) (rank-deficient noise). Rank-deficient noise is rejected at construction with a clear error message naming the limitation and the workarounds (ε-regularization or hand-rolled `FreidlinWentzellHamiltonian`); the error string does **not** cite the upstream issue number (PR/issue refs ages poorly and tie the library to GitHub state). The design that supports rank-deficient noise (`_integrate_kernel_drift!`, range/kernel projection in `update_p!`, range-only interpolation) is documented in this spec's Appendix A as a known follow-up so a future PR doesn't have to re-derive it from scratch.
**Replaces:** [PR #320](https://github.com/JuliaDynamics/CriticalTransitions.jl/pull/320) (draft).

## 1. Goals and non-goals

### Goals
* **(#328)** The auto-built `H_p` actually uses the SDE's covariance, not a hidden `Σ = I`. The trace-normalized convention from PR #329 is preserved.
* **(#275)** `fw_action`, `geometric_action`, `om_action`, and `minimize_geometric_action` all support state-dependent diffusion `σ(x)`. The case is enabled for both algorithms (gMAM and sgMAM).
* **(#326)** A single user-facing entry point `minimize_geometric_action(...)`. Dispatch is by input type:
  * `CoupledSDEs` ⇒ gMAM
  * `FreidlinWentzellHamiltonian` ⇒ sgMAM
* **(#263)** `minimize_geometric_action` no longer accepts `CoupledODEs`. The deterministic least-action path on an ODE is not a Freidlin-Wentzell rate object, and the existing `::ContinuousTimeDynamicalSystem` overload silently called `covariance_matrix(sys)` (which is undefined for `CoupledODEs`, the error reported in #263). Users wanting a drift-only minimum-action path use the `geometric_action(b::Function, path, arclength; A)` overload or the Hamiltonian path via `FreidlinWentzellHamiltonian{IIP, D}(H_x, H_p; a)`.
* Rename `ExtendedPhaseSpace` ⇒ `FreidlinWentzellHamiltonian`. No deprecation alias; the rename is a breaking change.
* Add a Hamiltonian-picture section to `docs/src/man/largedeviations.md` that motivates eq 17 of Grafke and Vanden-Eijnden 2018 and explains the three noise-shape branches.

### Non-goals
* Non-autonomous noise: still rejected.
* Non-Gaussian noise (jumps, Lévy): not supported; mentioned in docs only.
* AD-based `∂ₓa`: deferred. Use central finite differences.
* Rank-deficient noise (either constant or state-dependent): out of scope. The classifier detects it at construction and throws `ArgumentError` with a message naming the limitation and the workarounds (ε-regularization, or hand-rolled `FreidlinWentzellHamiltonian` for sgMAM). The error message does not cite issue or PR numbers.

## 2. Architecture

### 2.1 `FreidlinWentzellHamiltonian`

```julia
abstract type NoiseShape end
struct AdditiveNoise <: NoiseShape end
struct DiagonalNoise <: NoiseShape end
struct GeneralNoise  <: NoiseShape end

struct FreidlinWentzellHamiltonian{IIP, D, Hx, Hp, AF, NS <: NoiseShape}
    H_x::Hx          # (x, p) -> ∂H/∂x
    H_p::Hp          # (x, p) -> ∂H/∂p
    a::AF            # callable: x -> a(x), trace-normalized at u₀
end
```

Rank-deficient noise is detected at construction and rejected via an `ArgumentError` that names the limitation and the available workarounds. The supporting design is preserved in Appendix A.

`NS` is encoded as a **type parameter** (not a field), so `update_p!` and `update_x!` specialize by method dispatch with zero runtime branch cost.

### 2.2 Constructors

```julia
# auto from a CoupledSDEs / CoupledODEs
FreidlinWentzellHamiltonian(ds::ContinuousTimeDynamicalSystem)

# hand-rolled Hamiltonian
FreidlinWentzellHamiltonian{IIP, D}(
    H_x::Function, H_p::Function;
    a = Returns(LinearAlgebra.Diagonal(ones(D))),
)
```

The auto constructor:
1. Validates the SDE via `proper_FW_system(ds)` (autonomous noise only).
2. Pulls `σ(x) = diffusion_function(ds)`.
3. Builds `a(x) = σ(x) σ(x)ᵀ`.
4. Trace-normalizes at the reference state `u₀ = current_state(ds)` (Convention B): `s = tr(a(u₀))/D`; the stored `a` is `a(x)/s`. For an isotropic additive SDE with `noise_strength = σ_eff`, this gives the constant `I` and matches PR #329 exactly.
5. Classifies the noise shape at construction by sampling `a` at probe points around the reference state `u₀ = current_state(ds)` (or `zeros(D)` for the hand-rolled constructor): `u₀`, `u₀ + h·eᵢ`, `u₀ − h·eᵢ` for `i ∈ 1:D`, with `h = max(√eps(eltype(u₀)), 1e-6)`. If `cond(a(probe))` exceeds the rank tolerance at any probe, throw `ArgumentError` with a message naming the limitation and the available workarounds (ε-regularization, or hand-rolled `FreidlinWentzellHamiltonian` for sgMAM). Otherwise the resulting `NS` is:
   * **`AdditiveNoise`** if `a(probe_i) ≈ a(u₀)` to within `rtol = 1e-12` at every probe (no state dependence).
   * **`DiagonalNoise`** if every `a(probe_i)` is diagonal (off-diagonals within `atol = 1e-12`) and at least one probe differs from `a(u₀)`.
   * **`GeneralNoise`** if any `a(probe_i)` has non-zero off-diagonal entries.
6. Builds `H_p(x, p)` and `H_x(x, p)`:
   * Constant `a`: `H_p[:, i] = a*p[:, i] + b(x[:, i])` (broadcast over `Nt`).
   * Diagonal state-dependent: `H_p[:, i] = a_diag(x[:, i]) .* p[:, i] + b(x[:, i])`.
   * General: `H_p[:, i] = a(x[:, i]) * p[:, i] + b(x[:, i])`.
   * `H_x[:, i] = jac(x[:, i])ᵀ p[:, i] + ½ (∂ₓa(x[:, i]))·p[:, i]²` (FD for the `∂ₓa` term; zero for constant `a`).

The hand-rolled constructor trusts the user's `H_x`, `H_p`; it classifies `NS` by sampling the user's `a` at the same probe points.

### 2.3 Noise-shape classifier (shared between gMAM and sgMAM)

`CoupledSDEs` carries the noise structure through:

* `sys.noise_type :: NamedTuple{(:additive, :autonomous, :linear, :invertible), Bool}`, set at construction by `find_noise_type`. Note `:additive` is the *state-independent* flag (in DSB's naming, "additive" means state-independent); `:autonomous` is the time-independent flag. There is no `:state_independent` field.
* `diffusion_function(sys) :: (u, p, t) -> σ`, the callable.
* `covariance_matrix(sys)` returns `Σ = σσᵀ` only when `sys.noise_type[:invertible]` is **true**. For state-dependent noise (and for additive rank-deficient noise, which we reject) `covariance_matrix(sys) == nothing` (verified in DSB `core_systems/continuous_time_sde.jl`). The classifier therefore *never* calls `covariance_matrix(sys)` in the unsafe path; it builds `Σ` directly from `diffusion_function(sys)` so that the same code handles every input the classifier accepts.

These are runtime values, not type parameters. The classifier turns them into a singleton instance:

```julia
function _classify_noise_shape(sys::CoupledSDEs)
    sys.noise_type[:autonomous] ||
        throw(ArgumentError("non-autonomous noise not supported"))

    σ_fn = diffusion_function(sys)
    u₀   = current_state(sys)
    ps   = current_parameters(sys)

    a_of(u) = let σx = σ_fn(u, ps, 0.0)
        σ_mat = σx isa AbstractMatrix ? σx : LinearAlgebra.Diagonal(σx)
        σ_mat * σ_mat'
    end

    if sys.noise_type[:additive]
        # Constant a; no sampling needed. Build a directly from σ at u₀ rather than
        # going through `covariance_matrix(sys)`, which returns `nothing` when
        # `:invertible == false`.
        a = a_of(u₀)
        _is_rank_deficient(a) && throw(
            ArgumentError("rank-deficient noise is not supported. Workarounds: add a small ε on the noiseless variable to make the covariance invertible, or supply a Hamiltonian directly via FreidlinWentzellHamiltonian{IIP, D}(H_x, H_p)."),
        )
        return LinearAlgebra.isdiag(a) ? AdditiveNoise() : GeneralNoise()
    else
        # State-dependent: sample a(x) at probe points around u₀.
        probes = _probe_points(u₀)
        as = map(a_of, probes)

        any(_is_rank_deficient, as) && throw(
            ArgumentError("rank-deficient noise is not supported. Workarounds: add a small ε on the noiseless variable to make the covariance invertible, or supply a Hamiltonian directly via FreidlinWentzellHamiltonian{IIP, D}(H_x, H_p)."),
        )
        return all(LinearAlgebra.isdiag, as) ? DiagonalNoise() : GeneralNoise()
    end
end
```

Tolerances:

* `_is_rank_deficient(M)`: `cond(M) > 1/sqrt(eps(eltype(M)))` (or `det(M)` below `eps`).
* `_probe_points(u₀)`: returns `2D + 1` points `{u₀, u₀ ± h·eᵢ}` with `h = max(√eps(eltype(u₀)), 1e-6)`.

The same classifier is used by:

1. `FreidlinWentzellHamiltonian(ds)` constructor (Section 2.2): stores the result as a type parameter `NS`.
2. `minimize_geometric_action(::CoupledSDEs, ...)` entry point (Section 3.3): classifies, then constructs `GeometricGradientWorkspace{..., NS}` so subsequent `geometric_gradient_step!` calls dispatch on `NS` at compile time. After workspace construction there is no further runtime branching on the noise shape.

For the hand-rolled `FreidlinWentzellHamiltonian{IIP, D}(H_x, H_p; a)` constructor, the classifier is called on the user-supplied `a` callable (same `_probe_points` strategy around `zeros(D)`).

### 2.4 `Base.show` and `prettyprint`

`Base.show(io, sys::FreidlinWentzellHamiltonian)` displays the system dimension, in-place flag, and noise-shape tag (e.g. `"Freidlin-Wentzell Hamiltonian on 2-dimensional state space (DiagonalNoise) with out-of-place H_x and H_p"`). `prettyprint` is renamed accordingly and kept as an internal helper.

### 2.5 Removals (breaking)

The following names are removed without alias:

* `ExtendedPhaseSpace` (replaced by `FreidlinWentzellHamiltonian`).
* `minimize_simple_geometric_action` (replaced by `minimize_geometric_action(::FreidlinWentzellHamiltonian, ...)`).
* `proper_sgMAM_system` (the old sgMAM-specific check; replaced by `proper_FW_system`).
* `minimize_geometric_action(::CoupledODEs, ...)` (closes #263). The Hamiltonian constructor `FreidlinWentzellHamiltonian(::CoupledODEs)` is kept (a Hamiltonian with `a ≡ I` is well-defined for any drift) so users who want sgMAM on a deterministic system can still go via `minimize_geometric_action(FreidlinWentzellHamiltonian(ds), ...)`.

`proper_MAM_system` is kept (used by the OM/MAM path which still requires invertible additive noise).

## 3. Algorithm internals

Both algorithms (gMAM in `minimize_geometric_action.jl`, sgMAM in `sgMAM.jl`) get the same per-`NS` dispatch.

### 3.1 sgMAM `update_p!`

For all shapes, `b_ = H_p(x, 0)` (drift at every grid point).

#### `AdditiveNoise`
At construction we cache `a_diag` and `a_inv_diag`. Per call:
```
λ²(t) = Σᵢ(b_[i,t]²·a_inv_diag[i]) / Σᵢ(xdot[i,t]²·a_inv_diag[i])
p[i, t] = (λ(t)·xdot[i,t] − b_[i,t]) · a_inv_diag[i]
```
This collapses to the original `λ²(t) = ‖b‖²/‖xdot‖²`, `p = λ·xdot − b` when `a = I`.

#### `DiagonalNoise`
Per `t`: `a_t = a(x[:, t])`; `a_diag_t = diag(a_t)`. Same formula but `a_inv_diag` is per-`t`.

#### `GeneralNoise`
Per `t`: `a_t = a(x[:, t])`; LU-factor once, solve `a_t \ b_t` and `a_t \ xdot_t`.
```
λ²(t) = ⟨b_t, a_t⁻¹·b_t⟩ / ⟨xdot_t, a_t⁻¹·xdot_t⟩
p[:, t] = a_t⁻¹·(λ(t)·xdot_t − b_t)
```

_(`RankDeficientNoise` is rejected at construction; the design for handling it appears in Appendix A.)_

### 3.2 sgMAM `update_x!`

#### `AdditiveNoise`
Single shared `Tridiagonal(dl, d, du)` with `d = 1 + 2ϵλ²·c` and `c = a_inv_diag[dof]` constant. For isotropic `a = I` (the legacy fast path), `c = 1` and the code collapses exactly to today's implementation. Per-thread `LinearSolve` cache, reused across DOFs (current behavior).

#### `DiagonalNoise`
Precompute `a_at[Nx, Nt]` by evaluating `a(x[:, t])` once per timepoint at the start of the call. Per-DOF tridiagonal with `α[t] = ϵ·λ[t]²/a_at[dof, t]`. Threaded over DOFs. The cache strategy mirrors the additive path: one `LinearSolve` cache per thread, reused across DOFs via `reinit!` to swap in the new `Tridiagonal` and RHS. Concretely: at workspace construction we build `nthreads` caches with placeholder `Tridiagonal{T, Vector{T}}` of length `Nt-2`; each thread's inner loop over DOFs reuses its cache, overwriting `du`, `d`, `dl`, `b` in place before `solve!`.

#### `GeneralNoise`
Sparse block-tridiagonal of size `(Nt-2)·D × (Nt-2)·D`. Diagonal blocks `a(x[:, t]) + 2ϵλ(t)²·I`; off-diagonal blocks `-ϵλ(t)²·I`. The sparsity pattern is fixed across iterations:
* At construction, build the COO-to-`nzval` index map once and the placeholder `SparseMatrixCSC`.
* Per iteration, overwrite `nzval` in place using the cached map.
* Reuse a single `LinearSolve` cache so the symbolic factorization is shared across iterations.
* The cache is type-stable: the `SparseMatrixCSC{T, Int}` element type is fixed at construction (parameterized by `eltype(x_initial)`); the `LinearSolve.LinearCache` concrete type is computed once and stored in the workspace.

The pre-iteration work (evaluating `a(x[:, t])` for `t ∈ 1:Nt`, constructing the explicit RHS that includes the FD `∂ₓa·p²` term) is threaded over `t`. The linear solve itself is serial but uses BLAS internally; for `D ≤ 6` the per-iteration cost is dominated by the explicit RHS, not the solve.

This addresses PR #320's two performance TODOs and pulls "type-stable threaded LinearSolve caches" into scope.

### 3.3 gMAM (`minimize_geometric_action.jl`)

`GeometricGradientWorkspace` gains an `a_func::AF` field and an `NS<:NoiseShape` type parameter:

```julia
struct GeometricGradientWorkspace{Tupdate, ..., AF, NS<:NoiseShape}
    update::Tupdate
    a_func::AF
    # ... existing fields ...
end
```

`NS` is determined by calling `_classify_noise_shape(sys)` (Section 2.3) once when `minimize_geometric_action(::CoupledSDEs, ...)` constructs the workspace. After that, `geometric_gradient_step!` dispatches on `NS` via methods (compile-time specialization, no runtime branch per iteration):

* **`AdditiveNoise`**: existing fast path with constant `A = inv(a)`.
* **`DiagonalNoise`** / **`GeneralNoise`**: per-grid-point `a(x_i)` evaluation; explicit RHS uses `θ_i = a(x_i)⁻¹(λᵢ·φ'ᵢ − bᵢ)` and includes the `∂_x a` contribution to `H_φ` via central FD on `a`. The **implicit step stays a single shared tridiagonal** (independent of `a` per Heymann/Vanden-Eijnden 2008 eq 3.6 and Grafke et al. 2017 eq 19); state dependence enters only through the explicit RHS.

### 3.4 Performance commitments (in scope)

The implementation is held to these targets, all benchmarked as part of the PR:

* **Additive Maier-Stein (existing baseline)**: descent iteration cost within `1.05×` of the current `main` implementation. Acceptance threshold for the `benchmarks/kpo.jl` and `benchmarks/implementation_benchmarks/*.jl` suites.
* **Type-stable workspaces**: `GeometricGradientWorkspace` and `FreidlinWentzellHamiltonian` are fully type-stable end-to-end. Verified via `@code_warntype` on a hot path (`geometric_gradient_step!` / `update_x!`) showing no `Any` or `Union` returns. Added as a one-shot `Test.@inferred` assertion in the test suite for one path per branch.
* **Threaded LinearSolve caches** (additive + diagonal sgMAM `update_x!`): `nthreads` caches built at workspace construction, reused across DOFs via `LinearSolve.reinit!`. No allocations inside the inner DOF loop after warmup. Allocation-free path verified by `@allocations` regression test.
* **Sparse pattern cache (general sgMAM `update_x!`)**: COO-to-`nzval` index map built once at construction; `nzval` overwritten in place per iteration; one `LinearCache` reused. No reallocation of the sparse matrix in the inner iteration.
* **`a(x)` and `∂ₓa(x)` evaluation**: precomputed once per outer iteration (matrix `a_at[Nx, Nt]` or vector-of-matrices `a_at_full[Nt]`), then reused inside `update_x!` and `H_x`. Threaded over `t`. Avoids re-evaluating `a` at the same grid point inside the FD stencil and the implicit-step build.
* **Benchmark suite update**: add `benchmarks/multiplicative_noise.jl` exercising the three `NS` paths on the 1D OU multiplicative system and the 2D off-diagonal system; integrate into the existing `benchmark/` driver.

If any of these targets is missed by more than the stated threshold at PR time, the PR holds until either (a) the regression is understood and accepted with a written reason, or (b) the implementation is reworked.

### 3.5 Finite-difference `∂ₓa`

For `DiagonalNoise` and `GeneralNoise` (sgMAM `H_x` and gMAM explicit RHS), `∂ₓa` is computed by central FD:
* `h = max(√eps(eltype(x)), 1e-8)`.
* Per direction `l ∈ 1:D`: `∂_l a ≈ (a(x + h·eₗ) − a(x − h·eₗ)) / (2h)`.
* For `DiagonalNoise`, only the diagonal entries of `∂_l a` are used: `H_x[l] += ½ Σₖ pₖ²·(∂_l aₖₖ)`.
* For `GeneralNoise`, the full rank-3 tensor: `H_x[l] += ½ pᵀ(∂_l a)·p`.

The H_x closure built at construction encodes this choice based on `NS`.

## 4. API surface

### 4.1 Single entry point (#326)

```julia
# gMAM: tightened from ::ContinuousTimeDynamicalSystem to ::CoupledSDEs (closes #263).
# Internally classifies the noise shape of `sys`.
minimize_geometric_action(
    sys::CoupledSDEs,
    x_initial::Matrix,
    optimizer = GeometricGradient(; stepsize = 1.0);
    kwargs...,
)

# sgMAM: new dispatch, when the user supplies a Hamiltonian directly.
minimize_geometric_action(
    sys::FreidlinWentzellHamiltonian,
    x_initial::Matrix,
    optimizer = GeometricGradient(; stepsize = 1.0e3);
    kwargs...,
)
```

The existing `(sys, x_i, x_f, optimizer; npoints, ...)` overload that builds the path from endpoints is kept and is also tightened to `::CoupledSDEs`. `minimize_simple_geometric_action` is removed (no alias). Calling `minimize_geometric_action(::CoupledODEs, ...)` raises `MethodError`; the migration path is `geometric_action(drift, path, arclength; A)` for drift-only inference, or `FreidlinWentzellHamiltonian{IIP, D}(H_x, H_p)` for a hand-rolled Hamiltonian.

### 4.2 Action functions (#275)

`fw_action`, `geometric_action`, `om_action` get a small helper:

```julia
# Returns a callable x -> A(x), where A(x) is the inverse of the trace-normalized
# covariance a(x)/s with s = tr(a(u₀))/D. Rank-deficient noise is filtered out by
# _classify_noise_shape before we get here.
function _action_metric(sys::CoupledSDEs)
    _classify_noise_shape(sys)             # throws on rank-deficient or non-autonomous
    σ_fn = diffusion_function(sys)
    u₀ = current_state(sys); ps = current_parameters(sys)

    # `covariance_matrix(sys)` returns nothing whenever `:invertible == false`, so we
    # always build `a(x)` directly from `σ(x)` to keep the path uniform.
    a_of(x) = let σx = σ_fn(x, ps, 0.0)
        σ_mat = σx isa AbstractMatrix ? σx : LinearAlgebra.Diagonal(σx)
        σ_mat * σ_mat'
    end

    a₀ = a_of(u₀)
    s = tr(a₀) / size(a₀, 1)

    if sys.noise_type[:additive]
        A = inv(a₀ / s)
        return Returns(A)
    else
        return x -> inv(a_of(x) / s)
    end
end
```

The integrand evaluators (`fw_integrand`, `geometric_action`'s integrand) call the returned callable per grid point. For additive systems the callable is `Returns(A)`, so no per-point evaluation occurs and the existing fast path is preserved bit-for-bit.

### 4.3 `proper_*_system`

* New: `proper_FW_system(ds::CoupledSDEs)` requires only `:autonomous` noise (the rank-deficient rejection lives inside `_classify_noise_shape`). Used by `FreidlinWentzellHamiltonian(ds)` and by `minimize_geometric_action(::CoupledSDEs, ...)`.
* Existing: `proper_MAM_system` keeps the `(:additive, :invertible, :autonomous)` checks. Used by the finite-time `minimize_action` / OM path.
* Removed: `proper_sgMAM_system` (its `isdiag` check is replaced by the `NS` classifier).

### 4.4 Exports

`src/CriticalTransitions.jl` exports:

* New: `FreidlinWentzellHamiltonian`, `NoiseShape`, `AdditiveNoise`, `DiagonalNoise`, `GeneralNoise`.
* Removed: `ExtendedPhaseSpace`, `minimize_simple_geometric_action`.
* Existing: `minimize_geometric_action`, `fw_action`, `om_action`, `action`, `geometric_action`, `MinimumActionPath` (unchanged).

## 5. Testing strategy

### 5.1 Regression guards (must pass without modification)

* All existing cases in `test/largedeviations/sgMAM.jl` exercise `AdditiveNoise`. They must continue to pass when run via `FreidlinWentzellHamiltonian` and `minimize_geometric_action`.
* PR #329 conventions: `test/largedeviations/action_canonicalization.jl`.
* Maier-Stein β=1 analytic action `S = 1/2` via both gMAM and sgMAM.
* gMAM second-order convergence in `N`.

### 5.2 New: `test/largedeviations/multiplicative_noise.jl`

* **1D OU analytic baseline (DiagonalNoise):** `dX = −X dt + √(1+αX²) dW` on `[1, 0]`, `T = 1`. `fw_action` vs analytic Simpson; `rtol = 1e-4`.
* **Constant-`a` cross-check:** isotropic `noise_strength = σ` system vs `g = σ·I` system. Same fixed path; actions agree to `rtol = 1e-10`.
* **Additive limit:** 1D OU with `α = 1e-4`; gMAM and sgMAM both recover `S = 1` (additive limit) to `atol ≈ 3e-2`.
* **gMAM diagonal multiplicative vs Adam:** `GeometricGradient` and `Adam(0.01)` agree to `rtol = 0.05`.
* **gMAM general multiplicative vs Adam:** 2D off-diagonal `σ(x)`, `rtol = 0.1`.
* **Additive constant non-diagonal `Σ`:** Maier-Stein in a rotated basis with `covariance = R·diag(0.5, 2.0)·Rᵀ` for some rotation `R`. The classifier must return `GeneralNoise()` because the off-diagonal entries are non-zero. The action must agree (to `rtol = 1e-3`) with running the same system in the basis where `Σ` is diagonal (where it would classify as `DiagonalNoise`). Regression guard that the constant-non-diagonal path is exercised end-to-end and matches the diagonalized form.
* **sgMAM diagonal multiplicative vs Adam:** mirror via `FreidlinWentzellHamiltonian(ds)`; `rtol = 0.05`.
* **sgMAM general multiplicative vs Adam:** mirror; `rtol = 0.1`.
* **`NoiseShape` classification:** for each system, verify `FreidlinWentzellHamiltonian(ds)` ends up with the expected `NS` type parameter (introspect via `typeof`).

### 5.3 New: `test/largedeviations/rank_deficient_rejection.jl`

* **2D Langevin** `dx = p dt; dp = −∇V(x) dt − γp dt + √(2γT) dW`, `V(x) = (x²−1)²/4`. Constructing `FreidlinWentzellHamiltonian(ds)` throws `ArgumentError`. Calling `minimize_geometric_action(ds, ...)` (gMAM path) throws the same. Regression guard ensuring the rejection is informative.
* Error-message content: contains the strings `"rank-deficient"` and `"FreidlinWentzellHamiltonian"` (the latter cites the workaround). Prevents accidental degradation to a cryptic LinearAlgebra error and avoids tying the message to a GitHub issue number.

### 5.4 New: unified API tests (in `test/largedeviations/unified_api.jl`)

* `minimize_geometric_action(::CoupledSDEs, ...)` and `minimize_geometric_action(FreidlinWentzellHamiltonian(ds), ...)` converge to the same action (`rtol = 1e-2`) on Maier-Stein.
* Removed symbols: `@test !isdefined(CriticalTransitions, :minimize_simple_geometric_action)`, `@test !isdefined(CriticalTransitions, :ExtendedPhaseSpace)`.

### 5.5 Coverage matrix

| Noise type                                    | gMAM test                       | sgMAM test                                    |
|-----------------------------------------------|---------------------------------|-----------------------------------------------|
| Additive isotropic (`Σ = c·I`)                | existing Maier-Stein β=1        | existing sgMAM Maier-Stein                    |
| Additive constant non-diagonal `Σ`            | new: rotated Maier-Stein vs rotated additive baseline | new: same via Hamiltonian            |
| Diagonal multiplicative `a(x)::Diagonal`      | new: 1D OU + 2D vs Adam         | new: same systems via Hamiltonian             |
| General multiplicative `a(x)::Matrix`         | new: 2D off-diag vs Adam        | new: 2D off-diag vs Adam                      |
| Rank-deficient (any flavor)                   | rejection test (no descent run) | rejection test (no descent run)               |

## 6. File-level change list

* `src/largedeviations/sgMAM.jl`: rename struct, add `NS` parameter, refactor `update_p!`/`update_x!` into per-`NS` methods, drop `proper_sgMAM_system`.
* `src/largedeviations/minimize_geometric_action.jl`: add `a_func` and `NS` to workspace, refactor `geometric_gradient_step!` per `NS`, add `minimize_geometric_action(::FreidlinWentzellHamiltonian, ...)` wrapper.
* `src/largedeviations/action.jl`: add `_action_metric` helper; refactor `fw_integrand` and `geometric_action`'s integrand to use it; remove the hard `inv(Q)` precondition.
* `src/largedeviations/utils.jl`: add `proper_FW_system`; helpers for `NS` classification (`_classify_noise_shape`, `_is_diagonal_a`, `_constant_kernel_basis_or_throw`).
* `src/largedeviations/string_method.jl`: rename `ExtendedPhaseSpace` references to `FreidlinWentzellHamiltonian`. No semantic change (`string_method` calls `sys.H_p` at `p = 0`, which gives drift regardless of noise shape).
* `src/sde_utils.jl`: export `diffusion_function` from the DSB extension.
* `src/CriticalTransitions.jl`: update exports (add `FreidlinWentzellHamiltonian`, `NoiseShape`, four singletons; remove `ExtendedPhaseSpace`, `minimize_simple_geometric_action`).
* `test/largedeviations/sgMAM.jl`: update struct name and entry-point name in existing tests.
* `test/largedeviations/string_method.jl`: rename four `ExtendedPhaseSpace` references to `FreidlinWentzellHamiltonian`.
* New tests: `test/largedeviations/multiplicative_noise.jl`, `test/largedeviations/rank_deficient_rejection.jl`, `test/largedeviations/unified_api.jl`.
* `test/runtests.jl`: include the new test files.
* `docs/src/man/largedeviations.md`: add Hamiltonian-picture subsection; update method table; document the four `NoiseShape` cases and the Convention-B trace normalization.
* `examples/sgMAM_KPO.jl`, `examples/backtracking_KPO.jl`, `benchmarks/**`: rename struct name and entry-point.
* New benchmark: `benchmarks/multiplicative_noise.jl` exercising the three `NS` paths on the 1D OU multiplicative system and the 2D off-diagonal system; included in the existing `benchmark/` driver so regressions show up in CI alongside the additive baseline.

## 7. Risks and migration

### Breaking changes (user-visible)
* `ExtendedPhaseSpace` removed: users must rename to `FreidlinWentzellHamiltonian`.
* `minimize_simple_geometric_action` removed: users must call `minimize_geometric_action`.
* Action functions on previously-rejected SDEs now succeed (non-tightening): no user code that worked before stops working; some code that errored before now returns numerical values.

### Performance risks
* `update_x!` in `GeneralNoise` branch needs the cached sparse pattern + LinearSolve cache. Implemented as part of the initial PR (see Section 3.4); benchmarked against PR #320's uncached version. If the benchmark shows worse-than-baseline cost on the additive path (the common case), the PR holds.
* `H_x` FD adds O(D) `a` calls per grid point per iteration. Acceptable for typical sgMAM use (D ≤ 6); benchmarked on the 2D off-diagonal multiplicative system to confirm. Mitigation if needed: pre-evaluate `a` at `x ± h·eᵢ` once per `t` and reuse for both the FD stencil and the implicit-step build.
* `_classify_noise_shape` samples `a` at `2D + 1` points at construction. One-shot cost; negligible vs descent loop.
* Type-stability regression in any of the hot paths is treated as a release blocker; covered by the `@code_warntype` and `@allocations` checks listed in Section 3.4.

### Out of scope (follow-ups)
* **Rank-deficient noise (#325):** full design preserved in Appendix A. Once implemented, the classifier will return a new `RankDeficientNoise` singleton instead of erroring out, gMAM gets the range/kernel split, and sgMAM either supports it via the same machinery or directs users to gMAM.
* AD-based `∂ₓa` via `DifferentiationInterface.jl` (FD is fine for `D ≤ 6`; the AD swap is a drop-in once we hit a case where FD accuracy matters).

## 8. References

* Grafke and Vanden-Eijnden, *Numerical computation of rare events via large deviation theory*, Chaos 29 063118 (2019), [arXiv:1812.00681](https://arxiv.org/abs/1812.00681). Hamiltonian formulation: eqs 17, 18, 28, 33, 35.
* Grafke, Schäfer, Vanden-Eijnden, *Long term effects of small random perturbations* (2017). gMAM descent: eq 19.
* Heymann and Vanden-Eijnden, *The geometric minimum action method: a least action principle on the space of curves* (2008). Implicit operator structure: eq 3.6.
* Internal: PR #329 (action conventions), PR #320 (multiplicative-noise draft, replaced by this PR), issue #282 (the spurious `/2` bug, closed by PR #329).

## Appendix A. Deferred design: rank-deficient noise (#325)

This appendix documents the design that supports rank-deficient noise (constant kernel) so a future PR can pick it up without re-deriving. Reading this is not required for the present PR.

### A.1 Motivating cases
Second-order Langevin `dx = p dt; dp = (-∇V(x) - γp) dt + √(2γT) dW` has `a = diag(0, 2γT)`, rank `D-1` with kernel `span(e_x)`. Driven-dissipative quantum systems, fast/slow mode separation, and any model with "noise on a subset of variables" share the same structure: a constant linear subspace of the state space carries no noise, and on that subspace the path is forced to follow the deterministic drift exactly.

### A.2 Why the existing algorithms cannot handle it without modification
The FW action functional `S = (1/2) ∫ ||ẋ - b||²_{a⁻¹} dt` is degenerate in `ker(a)`: deviations in the noiseless direction contribute zero to the action. The optimizer therefore has no gradient signal to enforce the *physical* constraint `ẋ_⊥ = b_⊥(x)`, and the minimizer is non-unique. We need to add the constraint outside the action gradient.

Additionally, sgMAM's `update_p!` algebraically recovers `p = a⁻¹(λ ẋ - b)`. The kernel-direction component is undefined under a literal inverse; under pinv it's set to zero, which is consistent with the Hamiltonian `H_p` formula but doesn't tell us anything about `θ_⊥`. The kernel component of `θ` is in fact determined by the *other* Hamilton equation `θ̇ = -∇_x H` integrated along the path; sgMAM doesn't currently do that integration.

### A.3 Classifier change
Add a `RankDeficientNoise` singleton. In `_classify_noise_shape`, replace the `throw` with:
```julia
if _is_rank_deficient(a₀)
    # constant-kernel restriction: ker(a(probe_i)) must equal ker(a(u₀)) at every probe
    _validate_constant_kernel(as)
    return RankDeficientNoise()
end
```
`_validate_constant_kernel(as)` computes `ker(as[1])` once via SVD; for every other sample it verifies the same basis annihilates `as[i]` to within `atol = 1e-12`. Otherwise throws "state-dependent kernel direction not supported".

### A.4 Workspace and metric
* `_action_metric` returns `pinv(a/s)` instead of `inv(a/s)` for `RankDeficientNoise`. Kernel-direction residual is annihilated.
* Workspace caches the orthogonal decomposition once: range basis `U_r` (columns spanning range(a)), kernel basis `U_k`, and the range-restricted diagonal `D_r` (eigenvalues of `U_rᵀ a U_r`).

### A.5 gMAM `geometric_gradient_step!` (RankDeficientNoise)
* Compute `θ_i = pinv(a(x_i))·(λᵢ φ'ᵢ - bᵢ)` (kernel of θ is zero by pinv).
* Explicit RHS is built in the range subspace only; kernel coordinates of the RHS contribute zero.
* Implicit step is the same shared tridiagonal as additive (independent of `a`).
* **Constraint enforcement:** after the LinearSolve, project the updated path by integrating `x_⊥(t) = x_⊥(t-1) + Δt · b_⊥(x(t))` along the new grid. The integration runs over the path in arclength coordinates, using `λ(s)` to convert.

### A.6 sgMAM `update_p!` (RankDeficientNoise)
```
b_∥ = U_rᵀ·b_t;   b_⊥ = U_kᵀ·b_t
ẋ_∥ = U_rᵀ·ẋ_t
λ²(t) = ⟨b_∥, D_r⁻¹·b_∥⟩ / ⟨ẋ_∥, D_r⁻¹·ẋ_∥⟩
p_∥(t) = D_r⁻¹·(λ·ẋ_∥ - b_∥)
p_t    = U_r·p_∥(t)              # zero kernel component for p
```
The kernel component of `θ` (the conjugate momentum) is *not* free; it is determined by `θ̇_⊥ = -(∂_x H)_⊥` integrated along the path. The simplest implementation that still gives the correct instanton is to set `θ_⊥ = 0` (consistent with the leading-order minimizer) and rely on the kernel-direction constraint on `x` to do the work. For higher-order accuracy a future iteration can integrate `θ_⊥` separately.

### A.7 sgMAM `update_x!` (RankDeficientNoise)
Block-tridiagonal restricted to range-direction components; kernel-direction updates are deterministic:
```
x_⊥(t) = x_⊥(t-1) + Δs · b_⊥(t)
```

### A.8 Interaction with `interpolate_path!`
The existing `interpolate_path!` uses Euclidean arclength over all DOFs and would smear kernel-direction values during redistribution. A new helper `_interpolate_range_path!(path, α, s, U_r)` computes arclength using range-coord components only, redistributes range coords by linear interpolation, then re-derives kernel coords from the redistributed range path via `_integrate_kernel_drift!`. The rank-deficient outer-loop sequence becomes:

```
update_x!(...)                               # range coords change
_interpolate_range_path!(x, α, s, U_r)       # redistribute range coords
_integrate_kernel_drift!(x, U_r, U_k, λ, b)  # re-derive kernel coords
_sgmam_refresh!(xdot, p, lambda, x, H_p)     # recompute lambda, p
```

### A.9 Tests when the follow-up lands
* 2D Langevin: analytic FW action `S = ΔV/T` between minima; sgMAM and gMAM agree to `rtol = 0.05`.
* Constant-kernel detection failure: artificial SDE where `ker(a(x))` rotates triggers `ArgumentError`.
* `proper_FW_system` admits rank-deficient; `proper_MAM_system` still throws.
* Cross-check that the LDT instanton on Langevin agrees with the ε-regularized gMAM run as ε goes to zero.
