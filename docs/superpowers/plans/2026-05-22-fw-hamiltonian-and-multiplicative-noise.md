# Freidlin-Wentzell Hamiltonian, multiplicative noise, unified entry point: Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Land the design in [docs/superpowers/specs/2026-05-22-fw-hamiltonian-and-multiplicative-noise-design.md](../specs/2026-05-22-fw-hamiltonian-and-multiplicative-noise-design.md). Closes #263, #275, #326, #328. Rank-deficient noise (#325) is rejected at construction; the supporting design is preserved in the spec's Appendix A.

**Architecture:** Introduce a `NoiseShape` type hierarchy (`AdditiveNoise`, `DiagonalNoise`, `GeneralNoise`) encoded as a type parameter on a renamed `FreidlinWentzellHamiltonian` struct and on a new `NS`-parameterized `GeometricGradientWorkspace`. A single `_classify_noise_shape(::CoupledSDEs)` classifier runs once at construction; subsequent inner loops dispatch at compile time. Both algorithms (gMAM in `minimize_geometric_action.jl`, sgMAM in `sgMAM.jl`) get per-`NS` paths; action functions are generalized via a `_action_metric` helper.

**Tech Stack:** Julia 1.10+, `DynamicalSystemsBase`, `LinearSolve`, `SparseArrays`, `Test`, `BenchmarkTools` (for the benchmark file).

**Reference paths (read before starting):**
- Spec: `docs/superpowers/specs/2026-05-22-fw-hamiltonian-and-multiplicative-noise-design.md`
- Existing sgMAM impl: `src/largedeviations/sgMAM.jl`
- Existing gMAM impl: `src/largedeviations/minimize_geometric_action.jl`
- Existing action funcs: `src/largedeviations/action.jl`
- Existing string method: `src/largedeviations/string_method.jl`
- Existing sde utils: `src/sde_utils.jl`
- Existing utils: `src/largedeviations/utils.jl`
- Module top-level: `src/CriticalTransitions.jl`
- Tests live in `test/largedeviations/`; runner at `test/runtests.jl`.

**Working conventions used throughout this plan:**
- Use `git status` before each commit. Stage only the files you intended to change (avoid `git add -A`).
- Commit messages: short imperative subject (under 70 chars), one-line body if needed. No AI attribution trailer.
- **Never run `git commit` without explicit user approval.** "Proceed" or "looks good" do not count; wait for "commit" / "push" / equivalent.
- Run `julia --project=. -e 'using Pkg; Pkg.test()'` for the full suite. Note: `test_args=...` in `Pkg.test` is currently a **no-op** in this repo because `test/runtests.jl` does not process `ARGS`. The `Pkg.test(; test_args=[...])` invocations throughout this plan therefore run the entire suite. For tight inner-loop iteration on a single test file, invoke it directly:
  ```bash
  julia --project=. test/largedeviations/<file>.jl
  ```
  This bypasses the runner and runs just that file's `@testset` blocks. The "expected failure" steps below assume you read the targeted testset's output; other tests in the suite will continue to pass.
- Where a test asserts numerical agreement with an analytic value, use `@test isapprox(x, expected; rtol=...)` with tolerances from the spec.
- Do not introduce em-dashes or en-dashes in any file (repo policy via pre-commit hook).

---

## File Structure

**New files:**

- `src/largedeviations/noise_shape.jl`: `NoiseShape` types, `_classify_noise_shape`, helpers `_probe_points`, `_is_rank_deficient`.
- `test/largedeviations/noise_shape.jl`: classifier behavior tests.
- `test/largedeviations/multiplicative_noise.jl`: gMAM and sgMAM convergence and analytic baselines for diagonal and general multiplicative.
- `test/largedeviations/rank_deficient_rejection.jl`: rejection tests for 2D Langevin.
- `test/largedeviations/unified_api.jl`: entry-point dispatch tests, removed-symbol assertions, #263 regression.
- `benchmarks/multiplicative_noise.jl`: multiplicative noise benchmark cases.

**Modified files:**

- `src/CriticalTransitions.jl`: import `diffusion_function`; add new exports; remove old exports.
- `src/largedeviations/sgMAM.jl`: rename struct, add `NS` type parameter, refactor `update_p!`, `update_x!`, H_p/H_x closures.
- `src/largedeviations/minimize_geometric_action.jl`: add `a_func`/`NS` to workspace, refactor `geometric_gradient_step!`, tighten dispatch to `::CoupledSDEs`.
- `src/largedeviations/action.jl`: add `_action_metric`, generalize integrand evaluators.
- `src/largedeviations/utils.jl`: add `proper_FW_system`.
- `src/largedeviations/string_method.jl`: rename `ExtendedPhaseSpace` references.
- `test/largedeviations/sgMAM.jl`: rename struct in existing tests; update internals references.
- `test/largedeviations/string_method.jl`: rename `ExtendedPhaseSpace` references.
- `test/runtests.jl`: include new test files.
- `examples/sgMAM_KPO.jl`, `examples/backtracking_KPO.jl`: rename struct.
- `benchmarks/kpo.jl`, `benchmarks/implementation_benchmarks/KPO_sgMAM.jl`, `benchmarks/implementation_benchmarks/evaluate_drift.jl`, `benchmarks/implementation_benchmarks/string_method.jl`: rename struct.
- `docs/src/man/largedeviations.md`: add Hamiltonian-picture subsection, update method table, document `NoiseShape`.

---

## Phase 1: Foundation (NoiseShape and classifier)

### Task 1: Import `diffusion_function`

**Files:**
- Modify: `src/CriticalTransitions.jl` (the `using DynamicalSystemsBase: ...` block)

- [ ] **Step 1: Add `diffusion_function` to the import list**

In `src/CriticalTransitions.jl`, change the `using DynamicalSystemsBase:` block to include `diffusion_function`:

```julia
using DynamicalSystemsBase:
    DynamicalSystemsBase,
    CoupledSDEs,
    CoupledODEs,
    dynamic_rule,
    initial_state,
    current_state,
    set_state!,
    trajectory,
    jacobian,
    ContinuousTimeDynamicalSystem,
    initial_parameters,
    current_parameter,
    current_parameters,
    set_parameters!,
    initial_time,
    integrator,
    referenced_sciml_prob,
    covariance_matrix,
    diffusion_matrix,
    diffusion_function
```

- [ ] **Step 2: Verify the package still loads**

```bash
julia --project=. -e 'using CriticalTransitions; println(CriticalTransitions.diffusion_function)'
```
Expected: prints `DynamicalSystemsBase.StochasticSystemsBase.diffusion_function` (no `UndefVarError`).

- [ ] **Step 3: Ask user for approval, then commit**

Show the diff (`git diff src/CriticalTransitions.jl`). After explicit user approval ("commit"/"push"/"go ahead"), run:

```bash
git add src/CriticalTransitions.jl
git commit -m "import diffusion_function from DynamicalSystemsBase"
```

---

### Task 2: Add `NoiseShape` types

**Files:**
- Create: `src/largedeviations/noise_shape.jl`
- Modify: `src/CriticalTransitions.jl` (include the new file, export the names)
- Test: `test/largedeviations/noise_shape.jl`
- Modify: `test/runtests.jl` (include the new test file)

- [ ] **Step 1: Write the failing test**

Create `test/largedeviations/noise_shape.jl`:

```julia
using CriticalTransitions
using Test

@testset "NoiseShape hierarchy" begin
    @test NoiseShape === CriticalTransitions.NoiseShape
    @test AdditiveNoise <: NoiseShape
    @test DiagonalNoise <: NoiseShape
    @test GeneralNoise  <: NoiseShape
    @test AdditiveNoise() isa NoiseShape
    @test DiagonalNoise() isa NoiseShape
    @test GeneralNoise() isa NoiseShape
end
```

Wire it into the runner. In `test/runtests.jl`, add this line inside the `@testset "Large Deviations"`:

```julia
    include("largedeviations/noise_shape.jl")
```

- [ ] **Step 2: Run, expect failure**

```bash
julia --project=. -e 'using Pkg; Pkg.test(; test_args=["noise_shape"])'
```
Expected: `UndefVarError: NoiseShape not defined`.

- [ ] **Step 3: Implement**

Create `src/largedeviations/noise_shape.jl`:

```julia
"""
    NoiseShape

Abstract supertype for the diffusion-tensor shape tags used by the Freidlin-Wentzell
machinery. The tag is a singleton type carried as a type parameter so that the per-shape
algorithm branches dispatch at compile time.

Singletons:

* `AdditiveNoise`: `a(x)` is constant (state-independent) and `isdiag(a)` is true.
* `DiagonalNoise`: `a(x)` is diagonal at every probe point but varies with `x`.
* `GeneralNoise`: `a(x)` has off-diagonal entries (constant or state-dependent).

Rank-deficient noise (a singular diffusion tensor) is not supported in this release.
"""
abstract type NoiseShape end

struct AdditiveNoise <: NoiseShape end
struct DiagonalNoise <: NoiseShape end
struct GeneralNoise  <: NoiseShape end
```

In `src/CriticalTransitions.jl`, add `include("largedeviations/noise_shape.jl")` to the `include` block (find the existing `include("largedeviations/...")` lines and group the new file with them; the current includes are inside the module body after `include("sde_utils.jl")`). Then add to the `export` section:

```julia
export NoiseShape, AdditiveNoise, DiagonalNoise, GeneralNoise
```

- [ ] **Step 4: Run, expect pass**

```bash
julia --project=. -e 'using Pkg; Pkg.test(; test_args=["noise_shape"])'
```
Expected: all 6 `@test`s pass.

- [ ] **Step 5: Ask for approval, then commit**

```bash
git add src/largedeviations/noise_shape.jl src/CriticalTransitions.jl test/largedeviations/noise_shape.jl test/runtests.jl
git commit -m "add NoiseShape type hierarchy"
```

---

### Task 3: `_is_rank_deficient` and `_probe_points`

**Files:**
- Modify: `src/largedeviations/noise_shape.jl`
- Test: `test/largedeviations/noise_shape.jl` (add tests)

- [ ] **Step 1: Add failing tests**

Append to `test/largedeviations/noise_shape.jl`:

```julia
using LinearAlgebra

@testset "_is_rank_deficient" begin
    @test !CriticalTransitions._is_rank_deficient(LinearAlgebra.I(3))
    @test !CriticalTransitions._is_rank_deficient(LinearAlgebra.Diagonal([1.0, 2.0, 3.0]))
    @test  CriticalTransitions._is_rank_deficient(LinearAlgebra.Diagonal([1.0, 0.0, 3.0]))
    @test  CriticalTransitions._is_rank_deficient([1.0 1.0; 1.0 1.0])
end

@testset "_probe_points" begin
    u₀ = [0.5, -0.3]
    probes = CriticalTransitions._probe_points(u₀)
    @test length(probes) == 2 * length(u₀) + 1
    @test probes[1] == u₀
    h = max(sqrt(eps(eltype(u₀))), 1e-6)
    @test probes[2] ≈ u₀ + [h, 0]
    @test probes[3] ≈ u₀ - [h, 0]
    @test probes[4] ≈ u₀ + [0, h]
    @test probes[5] ≈ u₀ - [0, h]
end
```

- [ ] **Step 2: Run, expect failure**

```bash
julia --project=. -e 'using Pkg; Pkg.test(; test_args=["noise_shape"])'
```
Expected: `UndefVarError: _is_rank_deficient not defined`.

- [ ] **Step 3: Implement**

Append to `src/largedeviations/noise_shape.jl`:

```julia
"""
    _is_rank_deficient(M)

Return `true` if the matrix `M` is numerically singular. Uses `cond(M)` against
`1/sqrt(eps(eltype(M)))` as the tolerance.
"""
function _is_rank_deficient(M::AbstractMatrix)
    T = eltype(M)
    tol = 1 / sqrt(eps(real(T)))
    return LinearAlgebra.cond(Matrix(M)) > tol
end

"""
    _probe_points(u₀)

Return `2D + 1` probe points around `u₀` (where `D = length(u₀)`): the center point,
plus `u₀ ± h·eᵢ` for each canonical direction `eᵢ`, with `h = max(√eps, 1e-6)`. Used by
`_classify_noise_shape` to sample `a(x)` for state-dependence and rank-deficiency
detection.
"""
function _probe_points(u₀::AbstractVector)
    D = length(u₀)
    T = eltype(u₀)
    h = max(sqrt(eps(real(T))), T(1e-6))
    probes = Vector{Vector{T}}(undef, 2D + 1)
    probes[1] = collect(u₀)
    @inbounds for i in 1:D
        plus  = collect(u₀)
        minus = collect(u₀)
        plus[i]  += h
        minus[i] -= h
        probes[2i]     = plus
        probes[2i + 1] = minus
    end
    return probes
end
```

`LinearAlgebra` is in scope via the top-level `using LinearAlgebra: ...` in `src/CriticalTransitions.jl` which loads first.

- [ ] **Step 4: Run, expect pass**

```bash
julia --project=. -e 'using Pkg; Pkg.test(; test_args=["noise_shape"])'
```
Expected: all tests pass.

- [ ] **Step 5: Ask for approval, then commit**

```bash
git add src/largedeviations/noise_shape.jl test/largedeviations/noise_shape.jl
git commit -m "add _is_rank_deficient and _probe_points helpers"
```

---

### Task 4: `_classify_noise_shape(::CoupledSDEs)`

**Files:**
- Modify: `src/largedeviations/noise_shape.jl`
- Test: `test/largedeviations/noise_shape.jl` (add tests)

- [ ] **Step 1: Add failing tests**

Append to `test/largedeviations/noise_shape.jl`:

```julia
using StaticArrays

@testset "_classify_noise_shape: CoupledSDEs" begin
    f_lin(u, p, t) = SA[-u[1], -u[2]]
    ds_iso = CoupledSDEs(f_lin, SA[0.0, 0.0]; noise_strength = 1.0)
    @test CriticalTransitions._classify_noise_shape(ds_iso) isa AdditiveNoise

    R = [cos(0.3) -sin(0.3); sin(0.3) cos(0.3)]
    Q_user = R * Diagonal([0.5, 2.0]) * R'
    ds_nd = CoupledSDEs(f_lin, SA[0.0, 0.0]; covariance = Q_user)
    @test CriticalTransitions._classify_noise_shape(ds_nd) isa GeneralNoise

    g1d(u, p, t) = SA[sqrt(1 + 0.3 * u[1]^2);;]
    f1d(u, p, t) = SA[-u[1]]
    ds_diagmult = CoupledSDEs(
        f1d, SA[1.0]; g = g1d, noise_prototype = SMatrix{1, 1}(0.0),
    )
    @test CriticalTransitions._classify_noise_shape(ds_diagmult) isa DiagonalNoise

    f2(u, p, t) = SA[-u[1], -u[2]]
    function g2(u, p, t)
        s11 = 1 + 0.2 * u[1]; s22 = 1 + 0.2 * u[2]
        s12 = 0.3 * u[2];     s21 = 0.3 * u[1]
        return @SMatrix [s11 s12; s21 s22]
    end
    ds_gen = CoupledSDEs(
        f2, SA[1.0, 0.0]; g = g2, noise_prototype = SMatrix{2, 2}(zeros(2, 2)),
    )
    @test CriticalTransitions._classify_noise_shape(ds_gen) isa GeneralNoise
end

@testset "_classify_noise_shape: rank-deficient rejection" begin
    function langevin(u, p, t)
        x, p_ = u
        return SA[p_, -x - 0.1 * p_]
    end
    g_langevin(u, p, t) = SA[0.0 0.0; 0.0 sqrt(0.2)]
    ds_langevin = CoupledSDEs(
        langevin, SA[0.0, 0.0]; g = g_langevin,
        noise_prototype = SMatrix{2, 2}(zeros(2, 2)),
    )
    err = try
        CriticalTransitions._classify_noise_shape(ds_langevin)
        nothing
    catch e
        e
    end
    @test err isa ArgumentError
    @test occursin("rank-deficient", err.msg)
end
```

- [ ] **Step 2: Run, expect failure**

```bash
julia --project=. -e 'using Pkg; Pkg.test(; test_args=["noise_shape"])'
```
Expected: `UndefVarError: _classify_noise_shape not defined`.

- [ ] **Step 3: Implement**

Append to `src/largedeviations/noise_shape.jl`:

```julia
"""
    _classify_noise_shape(ds::CoupledSDEs)

Sample `a(x) = σ(x)σ(x)ᵀ` at `_probe_points(current_state(ds))` and return the
appropriate `NoiseShape` singleton. Throws `ArgumentError` if the noise is not
autonomous or if `a` is rank-deficient at any probe.

Uses `diffusion_function(ds)` directly because `covariance_matrix(ds)` returns
`nothing` whenever `:invertible == false` (including for rank-deficient additive
noise, which we want to detect cleanly here).
"""
function _classify_noise_shape(ds::CoupledSDEs)
    ds.noise_type[:autonomous] ||
        throw(ArgumentError("non-autonomous noise not supported"))

    σ_fn = diffusion_function(ds)
    u₀   = current_state(ds)
    ps   = current_parameters(ds)

    a_of(u) = let σx = σ_fn(u, ps, 0.0)
        σ_mat = σx isa AbstractMatrix ? σx : LinearAlgebra.Diagonal(σx)
        σ_mat * σ_mat'
    end

    if ds.noise_type[:additive]
        a = a_of(u₀)
        _is_rank_deficient(a) && throw(
            ArgumentError(
                "rank-deficient noise is not supported. Workarounds: add a small ε on the noiseless variable to make the covariance invertible, or supply a Hamiltonian directly via FreidlinWentzellHamiltonian{IIP, D}(H_x, H_p).",
            ),
        )
        return LinearAlgebra.isdiag(a) ? AdditiveNoise() : GeneralNoise()
    else
        probes = _probe_points(u₀)
        as = map(a_of, probes)
        any(_is_rank_deficient, as) && throw(
            ArgumentError(
                "rank-deficient noise is not supported. Workarounds: add a small ε on the noiseless variable to make the covariance invertible, or supply a Hamiltonian directly via FreidlinWentzellHamiltonian{IIP, D}(H_x, H_p).",
            ),
        )
        return all(LinearAlgebra.isdiag, as) ? DiagonalNoise() : GeneralNoise()
    end
end

# `CoupledODEs` has no noise; treat it as identity additive noise for the FW Hamiltonian.
_classify_noise_shape(::CoupledODEs) = AdditiveNoise()
```

- [ ] **Step 4: Run, expect pass**

```bash
julia --project=. -e 'using Pkg; Pkg.test(; test_args=["noise_shape"])'
```
Expected: all tests pass.

- [ ] **Step 5: Ask for approval, then commit**

```bash
git add src/largedeviations/noise_shape.jl test/largedeviations/noise_shape.jl
git commit -m "add _classify_noise_shape with rank-deficient rejection"
```

---

### Task 5: `proper_FW_system`

**Files:**
- Modify: `src/largedeviations/utils.jl` (add new function)
- Test: `test/largedeviations/noise_shape.jl` (add a test)

- [ ] **Step 1: Add failing test**

Append to `test/largedeviations/noise_shape.jl`:

```julia
@testset "proper_FW_system" begin
    f_lin(u, p, t) = SA[-u[1], -u[2]]
    ds_iso = CoupledSDEs(f_lin, SA[0.0, 0.0]; noise_strength = 1.0)
    @test CriticalTransitions.proper_FW_system(ds_iso) === nothing

    g1d(u, p, t) = SA[sqrt(1 + 0.3 * u[1]^2);;]
    f1d(u, p, t) = SA[-u[1]]
    ds_mult = CoupledSDEs(
        f1d, SA[1.0]; g = g1d, noise_prototype = SMatrix{1, 1}(0.0),
    )
    @test CriticalTransitions.proper_FW_system(ds_mult) === nothing
end
```

- [ ] **Step 2: Run, expect failure**

```bash
julia --project=. -e 'using Pkg; Pkg.test(; test_args=["noise_shape"])'
```
Expected: `UndefVarError: proper_FW_system not defined`.

- [ ] **Step 3: Implement**

In `src/largedeviations/utils.jl`, add immediately after `proper_MAM_system`:

```julia
"""
    proper_FW_system(ds::CoupledSDEs)

Validates that `ds` is a valid input for the Freidlin-Wentzell Hamiltonian path
(gMAM or sgMAM via `FreidlinWentzellHamiltonian(ds)`). Only requires autonomous noise;
rank-deficiency is detected by `_classify_noise_shape` (called at workspace
construction). Returns `nothing` on success; throws `ArgumentError` otherwise.
"""
function proper_FW_system(ds::CoupledSDEs)
    if !ds.noise_type[:autonomous]
        throw(
            ArgumentError(
                "Freidlin-Wentzell methods are only applicable for autonomous noise.",
            ),
        )
    end
    return nothing
end
```

- [ ] **Step 4: Run, expect pass**

```bash
julia --project=. -e 'using Pkg; Pkg.test(; test_args=["noise_shape"])'
```
Expected: all tests pass.

- [ ] **Step 5: Ask for approval, then commit**

```bash
git add src/largedeviations/utils.jl test/largedeviations/noise_shape.jl
git commit -m "add proper_FW_system (autonomous-only check)"
```

---

## Phase 2: Rename struct and restructure

### Task 6: Rename `ExtendedPhaseSpace` to `FreidlinWentzellHamiltonian`

Mechanical multi-file rename. Existing semantics preserved.

**Files:**
- Modify: `src/largedeviations/sgMAM.jl`
- Modify: `src/largedeviations/string_method.jl`
- Modify: `src/CriticalTransitions.jl` (export rename)
- Modify: `test/largedeviations/sgMAM.jl`
- Modify: `test/largedeviations/string_method.jl`
- Modify: `examples/sgMAM_KPO.jl`, `examples/backtracking_KPO.jl`
- Modify: `benchmarks/kpo.jl`, `benchmarks/implementation_benchmarks/KPO_sgMAM.jl`, `benchmarks/implementation_benchmarks/evaluate_drift.jl`, `benchmarks/implementation_benchmarks/string_method.jl`
- Modify: `docs/src/man/largedeviations.md`, `docs/src/examples/sgMAM_KPO.md` (if it exists)

- [ ] **Step 1: Baseline**

```bash
julia --project=. -e 'using Pkg; Pkg.test(; test_args=["sgMAM", "string_method"])'
```
Expected: PASS.

- [ ] **Step 2: Enumerate occurrences**

```bash
git grep -nE 'ExtendedPhaseSpace' -- ':!docs/superpowers/'
```
Replace `ExtendedPhaseSpace` with `FreidlinWentzellHamiltonian` everywhere it appears in those files (struct declaration, type ascriptions, isa-checks, docstrings, `prettyprint`, and the corresponding `Base.show` method).

The export line in `src/CriticalTransitions.jl` should change from `export minimize_simple_geometric_action, ExtendedPhaseSpace` to `export minimize_simple_geometric_action, FreidlinWentzellHamiltonian` (the `minimize_simple_geometric_action` export is removed later in Task 23).

- [ ] **Step 3: Re-run tests**

```bash
julia --project=. -e 'using Pkg; Pkg.test(; test_args=["sgMAM", "string_method"])'
```
Expected: PASS (semantics unchanged).

- [ ] **Step 4: Ask for approval, then commit**

```bash
git status
git add -u
git commit -m "rename ExtendedPhaseSpace to FreidlinWentzellHamiltonian"
```

---

### Task 7: Add `a::AF` field and `NS` type parameter

**Files:**
- Modify: `src/largedeviations/sgMAM.jl`
- Modify: `src/largedeviations/noise_shape.jl` (add `_classify_user_a`)
- Test: `test/largedeviations/sgMAM.jl` (add type-parameter introspection test)

- [ ] **Step 1: Add failing test**

In `test/largedeviations/sgMAM.jl`, add:

```julia
@testset "FreidlinWentzellHamiltonian carries NoiseShape" begin
    f_lin(u, p, t) = SA[-u[1], -u[2]]
    ds_ode = CoupledODEs(f_lin, SA[0.0, 0.0])
    sys_ode = FreidlinWentzellHamiltonian(ds_ode)
    @test sys_ode isa FreidlinWentzellHamiltonian{<:Any, 2, <:Any, <:Any, <:Any, AdditiveNoise}
    @test sys_ode.a(zeros(2)) ≈ LinearAlgebra.Diagonal(ones(2))

    ds_iso = CoupledSDEs(f_lin, SA[0.0, 0.0]; noise_strength = 1.0)
    sys_iso = FreidlinWentzellHamiltonian(ds_iso)
    @test sys_iso isa FreidlinWentzellHamiltonian{<:Any, 2, <:Any, <:Any, <:Any, AdditiveNoise}

    H_x_user(x, p) = zeros(size(x))
    H_p_user(x, p) = ones(size(x))
    sys_user = FreidlinWentzellHamiltonian{false, 2}(H_x_user, H_p_user)
    @test sys_user isa FreidlinWentzellHamiltonian{false, 2, <:Any, <:Any, <:Any, AdditiveNoise}
end
```

- [ ] **Step 2: Run, expect failure**

```bash
julia --project=. -e 'using Pkg; Pkg.test(; test_args=["sgMAM"])'
```
Expected: existing tests still pass; new test fails because the struct has neither `a` field nor a sixth type parameter.

- [ ] **Step 3: Restructure the struct**

In `src/largedeviations/sgMAM.jl`, replace the struct definition and both constructors:

```julia
"""
A structure representing an extended phase space system where your dissipative vector
field is embedded in a doubled dimensional phase space. The struct stores the partial
derivatives `H_x`, `H_p` of the Freidlin/Wentzell Hamiltonian

```math
H(x, p) = ⟨b(x), p⟩ + (1/2) ⟨p, a(x)·p⟩,
```

a callable `a` for the diffusion tensor, and a `NoiseShape` type parameter `NS` that
encodes how the inner algorithm loops dispatch.
"""
struct FreidlinWentzellHamiltonian{IIP, D, Hx, Hp, AF, NS <: NoiseShape}
    H_x::Hx
    H_p::Hp
    a::AF
end

function FreidlinWentzellHamiltonian(ds::ContinuousTimeDynamicalSystem)
    D = dimension(ds)
    if ds isa CoupledSDEs
        proper_FW_system(ds)
    end
    NS = _classify_noise_shape(ds)

    σ_fn = ds isa CoupledSDEs ? diffusion_function(ds) : nothing
    ps   = current_parameters(ds)

    # Decision for `a(x)`:
    #  * No noise (CoupledODEs):                 a ≡ I, wrap in Returns.
    #  * Additive (state-independent) SDE:       a is constant; wrap the (trace-
    #    normalized) value in Returns regardless of diagonal-ness so that update_p!
    #    /update_x! do not repeatedly re-evaluate σ_fn for a constant tensor.
    #  * State-dependent SDE:                    closure that evaluates σ(x)σ(x)ᵀ
    #    per call, trace-normalized at u₀.
    is_constant_a = (ds isa CoupledODEs) || (ds isa CoupledSDEs && ds.noise_type[:additive])

    a_callable = if ds isa CoupledODEs
        Returns(LinearAlgebra.Diagonal(ones(Float64, D)))
    elseif is_constant_a
        u₀ = current_state(ds)
        σ0 = σ_fn(u₀, ps, 0.0)
        σ_mat = σ0 isa AbstractMatrix ? σ0 : LinearAlgebra.Diagonal(σ0)
        a0 = σ_mat * σ_mat'
        s = LinearAlgebra.tr(a0) / D
        # Preserve Diagonal type for diagonal Σ; this also keeps the AdditiveNoise
        # dispatch on update_x! fully type-stable.
        a_const = LinearAlgebra.isdiag(a0) ?
            LinearAlgebra.Diagonal(collect(LinearAlgebra.diag(a0)) ./ s) :
            Matrix(a0 ./ s)
        Returns(a_const)
    else
        u₀ = current_state(ds)
        σ0 = σ_fn(u₀, ps, 0.0)
        σ_mat0 = σ0 isa AbstractMatrix ? σ0 : LinearAlgebra.Diagonal(σ0)
        a0 = σ_mat0 * σ_mat0'
        s = LinearAlgebra.tr(a0) / D
        let σ_fn = σ_fn, ps = ps, s = s
            x -> begin
                σx = σ_fn(x, ps, 0.0)
                σ_mat = σx isa AbstractMatrix ? σx : LinearAlgebra.Diagonal(σx)
                (σ_mat * σ_mat') / s
            end
        end
    end

    f   = dynamic_rule(ds)
    jac = jacobian(ds)

    # Decide once whether the FD ∂ₓa term is needed. Constant a (additive SDE or
    # CoupledODEs) gives ∂ₓa ≡ 0, so we skip the FD section entirely; this also
    # handles the "constant non-diagonal Σ" case correctly even though it classifies
    # as GeneralNoise.
    skip_ax_fd = is_constant_a

    function H_x(x, p)
        Hx = similar(x)
        Nt = size(x, 2); Dn = size(x, 1)
        h_fd = max(sqrt(eps(eltype(x))), eltype(x)(1e-8))
        for idx in 1:Nt
            xi = x[:, idx]
            pi_v = p[:, idx]
            jax = jac(xi, ps, 0.0)
            @inbounds for idc in 1:Dn
                Hx[idc, idx] = dot(jax[:, idc], pi_v)
            end
            if !skip_ax_fd
                e = zeros(eltype(xi), Dn)
                sample = a_callable(xi)
                is_diag = sample isa LinearAlgebra.Diagonal || LinearAlgebra.isdiag(sample)
                @inbounds for l in 1:Dn
                    fill!(e, 0); e[l] = h_fd
                    ap = a_callable(xi .+ e)
                    am = a_callable(xi .- e)
                    if is_diag
                        dla = (LinearAlgebra.diag(ap) .- LinearAlgebra.diag(am)) ./ (2 * h_fd)
                        Hx[l, idx] += 0.5 * dot(pi_v .^ 2, dla)
                    else
                        Hx[l, idx] += 0.5 * dot(pi_v, ((ap .- am) ./ (2 * h_fd)) * pi_v)
                    end
                end
            end
        end
        return Hx
    end

    function H_p(x, p)
        Hp = similar(x)
        for idx in 1:size(x, 2)
            a_x = a_callable(x[:, idx])
            Hp[:, idx] = a_x * p[:, idx] .+ f(x[:, idx], ps, 0.0)
        end
        return Hp
    end

    return FreidlinWentzellHamiltonian{
        isinplace(ds), D, typeof(H_x), typeof(H_p), typeof(a_callable), typeof(NS),
    }(H_x, H_p, a_callable)
end

function FreidlinWentzellHamiltonian{IIP, D}(
        H_x::Function, H_p::Function;
        a = Returns(LinearAlgebra.Diagonal(ones(Float64, D))),
    ) where {IIP, D}
    NS = _classify_user_a(a, D)
    return FreidlinWentzellHamiltonian{IIP, D, typeof(H_x), typeof(H_p), typeof(a), typeof(NS)}(
        H_x, H_p, a,
    )
end
```

Note that the `H_x` closure already includes the `∂ₓa` finite-difference term (so Task 17 collapses into this task). Constant-`a` paths skip the FD section via the `NS isa AdditiveNoise` check.

Add the `_classify_user_a` helper to `src/largedeviations/noise_shape.jl`:

```julia
"""
    _classify_user_a(a_callable, D)

Classify a user-supplied `a(x)` callable. Probes `a` around `zeros(D)`; detects
constant-vs-state-dependent structurally (the callable is a `Base.Returns`) rather
than by numerical comparison, which avoids false positives on barely-state-dependent
functions.
"""
function _classify_user_a(a_callable, D::Int)
    probes = _probe_points(zeros(Float64, D))
    as = map(a_callable, probes)
    if any(_is_rank_deficient, as)
        throw(
            ArgumentError(
                "rank-deficient noise is not supported. Workarounds: add a small ε on the noiseless variable to make the covariance invertible, or supply a Hamiltonian directly via FreidlinWentzellHamiltonian{IIP, D}(H_x, H_p).",
            ),
        )
    end
    if a_callable isa Base.Returns
        return LinearAlgebra.isdiag(as[1]) ? AdditiveNoise() : GeneralNoise()
    end
    return all(LinearAlgebra.isdiag, as) ? DiagonalNoise() : GeneralNoise()
end
```

Update `prettyprint`:

```julia
function prettyprint(mlp::FreidlinWentzellHamiltonian{IIP, D, Hx, Hp, AF, NS}) where {IIP, D, Hx, Hp, AF, NS}
    return "Freidlin-Wentzell Hamiltonian on $D-dimensional state space ($(NS)) with $(IIP ? "in-place" : "out-of-place") H_x and H_p"
end
```

- [ ] **Step 4: Run tests**

```bash
julia --project=. -e 'using Pkg; Pkg.test(; test_args=["sgMAM", "string_method"])'
```
Expected: all existing sgMAM and string-method tests pass, plus the new NoiseShape test.

- [ ] **Step 5: Ask for approval, then commit**

```bash
git add src/largedeviations/sgMAM.jl src/largedeviations/noise_shape.jl test/largedeviations/sgMAM.jl
git commit -m "add a callable, NS type parameter, and FD H_x to FreidlinWentzellHamiltonian"
```

---

### Task 8: Remove `proper_sgMAM_system`

**Files:**
- Modify: `src/largedeviations/sgMAM.jl`

- [ ] **Step 1: Verify no other callers**

```bash
git grep -n 'proper_sgMAM_system' -- ':!docs/superpowers/'
```
Expected: only definition site in `src/largedeviations/sgMAM.jl` (after Task 7 it is unreferenced).

- [ ] **Step 2: Delete the function from `sgMAM.jl`**

- [ ] **Step 3: Run tests**

```bash
julia --project=. -e 'using Pkg; Pkg.test(; test_args=["sgMAM"])'
```
Expected: PASS.

- [ ] **Step 4: Ask for approval, then commit**

```bash
git add src/largedeviations/sgMAM.jl
git commit -m "remove proper_sgMAM_system (replaced by proper_FW_system)"
```

---

## Phase 3: sgMAM `update_p!` per `NS`

### Task 9: Refactor `update_p!` to dispatch on NS (AdditiveNoise preserves current behavior)

**Files:**
- Modify: `src/largedeviations/sgMAM.jl`
- Modify: `test/largedeviations/sgMAM.jl` (signature change)

- [ ] **Step 1: Run existing tests as baseline**

```bash
julia --project=. -e 'using Pkg; Pkg.test(; test_args=["sgMAM"])'
```
Expected: PASS.

- [ ] **Step 2: Refactor**

Replace the current `update_p!` body with a dispatching wrapper plus per-NS methods:

```julia
function update_p!(p, lambda, x, xdot, sys::FreidlinWentzellHamiltonian)
    _update_p!(p, lambda, x, xdot, sys)
    return nothing
end

function _update_p!(
        p, lambda, x, xdot,
        sys::FreidlinWentzellHamiltonian{IIP, D, Hx, Hp, AF, AdditiveNoise},
    ) where {IIP, D, Hx, Hp, AF}
    b_ = sys.H_p(x, zero(x))
    a_diag = LinearAlgebra.diag(sys.a(view(x, :, 1)))
    a_inv  = inv.(a_diag)
    num = sum(b_ .^ 2 .* a_inv; dims = 1)
    den = sum(xdot .^ 2 .* a_inv; dims = 1)
    lambda .= sqrt.(num ./ den)
    lambda[1] = 0
    lambda[end] = 0
    p .= (lambda .* xdot .- b_) .* a_inv
    return nothing
end
```

Update `_sgmam_refresh!`:

```julia
function _sgmam_refresh!(xdot, p, lambda, x, sys)
    central_diff!(xdot, x)
    update_p!(p, lambda, x, xdot, sys)
    return nothing
end
```

Update `update!`:

```julia
function update!(x, xdot, xdotdot, p, pdot, lambda, sys, ϵ)
    central_diff!(pdot, p)
    Hx = sys.H_x(x, p)
    central_diff!(xdotdot, xdot)
    return update_x!(x, lambda, pdot, xdotdot, Hx, sys, ϵ)
end
```

Update every call site in `minimize_simple_geometric_action` (both the GeometricGradient and AdaptiveGeometricGradient overloads). Replace:

- `_sgmam_refresh!(xdot, p, lambda, x, H_p)` with `_sgmam_refresh!(xdot, p, lambda, x, sys)`.
- `update!(x, xdot, xdotdot, p, pdot, lambda, H_x, H_p, ϵ)` with `update!(x, xdot, xdotdot, p, pdot, lambda, sys, ϵ)`.

Update the existing test in `test/largedeviations/sgMAM.jl` that calls `CT.update_p!(p, λ, x0, xdot, sys.H_p)` to `CT.update_p!(p, λ, x0, xdot, sys)`.

- [ ] **Step 3: Run, expect pass**

```bash
julia --project=. -e 'using Pkg; Pkg.test(; test_args=["sgMAM"])'
```
Expected: PASS.

- [ ] **Step 4: Ask for approval, then commit**

```bash
git add src/largedeviations/sgMAM.jl test/largedeviations/sgMAM.jl
git commit -m "refactor update_p! to NS-dispatched methods (additive preserved)"
```

---

### Task 10: `update_p!` for `DiagonalNoise`

**Files:**
- Modify: `src/largedeviations/sgMAM.jl`
- Create: `test/largedeviations/multiplicative_noise.jl`
- Modify: `test/runtests.jl`

- [ ] **Step 1: Create test file with a failing test**

Create `test/largedeviations/multiplicative_noise.jl`:

```julia
using CriticalTransitions, StaticArrays
using Test
using LinearAlgebra
using Random

const CT = CriticalTransitions

@testset "DiagonalNoise update_p!: 1D OU multiplicative" begin
    α = 0.3
    b1d(u, p, t) = SA[-u[1]]
    g1d(u, p, t) = SA[sqrt(1 + α * u[1]^2);;]
    ds = CoupledSDEs(b1d, SA[1.0]; g = g1d, noise_prototype = SMatrix{1, 1}(0.0))

    sys = FreidlinWentzellHamiltonian(ds)
    @test sys isa FreidlinWentzellHamiltonian{<:Any, 1, <:Any, <:Any, <:Any, DiagonalNoise}

    Nt = 40
    x0 = reshape(collect(range(1.0, -1.0; length = Nt)), 1, Nt)
    s_arc = range(0; stop = 1, length = Nt)
    α_arc = zeros(Nt)
    CT.interpolate_path!(x0, α_arc, s_arc)
    xdot = zeros(size(x0))
    p = zeros(size(x0))
    λ = zeros(1, Nt)
    CT.central_diff!(xdot, x0)
    CT.update_p!(p, λ, x0, xdot, sys)
    @test all(isfinite, λ)
    @test all(isfinite, p)
end
```

Add `include("largedeviations/multiplicative_noise.jl")` to `test/runtests.jl`.

- [ ] **Step 2: Run, expect failure**

```bash
julia --project=. -e 'using Pkg; Pkg.test(; test_args=["multiplicative_noise"])'
```
Expected: `MethodError` on the `DiagonalNoise` dispatch.

- [ ] **Step 3: Implement**

Append to `src/largedeviations/sgMAM.jl`:

```julia
function _update_p!(
        p, lambda, x, xdot,
        sys::FreidlinWentzellHamiltonian{IIP, D, Hx, Hp, AF, DiagonalNoise},
    ) where {IIP, D, Hx, Hp, AF}
    b_ = sys.H_p(x, zero(x))
    Nt = size(x, 2)
    for t in 1:Nt
        a_t   = sys.a(view(x, :, t))
        diag_t = LinearAlgebra.diag(a_t)
        inv_t  = inv.(diag_t)
        num_t = sum(b_[:, t] .^ 2 .* inv_t)
        den_t = sum(xdot[:, t] .^ 2 .* inv_t)
        λ_t = den_t > 1e-28 ? sqrt(num_t / den_t) : zero(eltype(b_))
        lambda[1, t] = isfinite(λ_t) ? λ_t : zero(eltype(b_))
        @inbounds for i in 1:D
            p[i, t] = (lambda[1, t] * xdot[i, t] - b_[i, t]) * inv_t[i]
        end
    end
    lambda[1, 1] = 0
    lambda[1, end] = 0
    return nothing
end
```

- [ ] **Step 4: Run, expect pass**

```bash
julia --project=. -e 'using Pkg; Pkg.test(; test_args=["multiplicative_noise"])'
```
Expected: PASS.

- [ ] **Step 5: Ask for approval, then commit**

```bash
git add src/largedeviations/sgMAM.jl test/largedeviations/multiplicative_noise.jl test/runtests.jl
git commit -m "add DiagonalNoise update_p! for state-dependent diagonal a(x)"
```

---

### Task 11: `update_p!` for `GeneralNoise`

**Files:**
- Modify: `src/largedeviations/sgMAM.jl`
- Modify: `test/largedeviations/multiplicative_noise.jl`

- [ ] **Step 1: Add failing test**

Append to `test/largedeviations/multiplicative_noise.jl`:

```julia
@testset "GeneralNoise update_p!: 2D off-diagonal" begin
    b2(u, p, t) = SA[-u[1], -u[2]]
    function g2(u, p, t)
        s11 = 1 + 0.2 * u[1]; s22 = 1 + 0.2 * u[2]
        s12 = 0.3 * u[2];     s21 = 0.3 * u[1]
        return @SMatrix [s11 s12; s21 s22]
    end
    ds = CoupledSDEs(
        b2, SA[1.0, 0.0]; g = g2, noise_prototype = SMatrix{2, 2}(zeros(2, 2)),
    )
    sys = FreidlinWentzellHamiltonian(ds)
    @test sys isa FreidlinWentzellHamiltonian{<:Any, 2, <:Any, <:Any, <:Any, GeneralNoise}

    Nt = 40
    xx = collect(range(1.0, 0.0; length = Nt))
    yy = collect(range(0.0, 1.0; length = Nt))
    x0 = Matrix([xx yy]')
    s_arc = range(0; stop = 1, length = Nt)
    α_arc = zeros(Nt)
    CT.interpolate_path!(x0, α_arc, s_arc)
    xdot = zeros(size(x0))
    p = zeros(size(x0))
    λ = zeros(1, Nt)
    CT.central_diff!(xdot, x0)
    CT.update_p!(p, λ, x0, xdot, sys)
    @test all(isfinite, λ)
    @test all(isfinite, p)
end
```

- [ ] **Step 2: Run, expect failure**

```bash
julia --project=. -e 'using Pkg; Pkg.test(; test_args=["multiplicative_noise"])'
```
Expected: `MethodError` on the `GeneralNoise` dispatch.

- [ ] **Step 3: Implement**

Append to `src/largedeviations/sgMAM.jl`:

```julia
function _update_p!(
        p, lambda, x, xdot,
        sys::FreidlinWentzellHamiltonian{IIP, D, Hx, Hp, AF, GeneralNoise},
    ) where {IIP, D, Hx, Hp, AF}
    b_ = sys.H_p(x, zero(x))
    Nt = size(x, 2)
    inv_b    = Vector{eltype(b_)}(undef, D)
    inv_xdot = Vector{eltype(b_)}(undef, D)
    rhs      = Vector{eltype(b_)}(undef, D)
    for t in 1:Nt
        a_t = sys.a(view(x, :, t))
        inv_b    .= a_t \ Vector(view(b_, :, t))
        inv_xdot .= a_t \ Vector(view(xdot, :, t))
        num_t = dot(view(b_, :, t),   inv_b)
        den_t = dot(view(xdot, :, t), inv_xdot)
        λ_t = den_t > 1e-28 ? sqrt(num_t / den_t) : zero(eltype(b_))
        lambda[1, t] = isfinite(λ_t) ? λ_t : zero(eltype(b_))
        @inbounds for i in 1:D
            rhs[i] = lambda[1, t] * xdot[i, t] - b_[i, t]
        end
        p[:, t] .= a_t \ rhs
    end
    lambda[1, 1] = 0
    lambda[1, end] = 0
    return nothing
end
```

- [ ] **Step 4: Run, expect pass**

```bash
julia --project=. -e 'using Pkg; Pkg.test(; test_args=["multiplicative_noise"])'
```
Expected: PASS.

- [ ] **Step 5: Ask for approval, then commit**

```bash
git add src/largedeviations/sgMAM.jl test/largedeviations/multiplicative_noise.jl
git commit -m "add GeneralNoise update_p! for state-dependent matrix a(x)"
```

---

## Phase 4: sgMAM `update_x!` per `NS`

### Task 12: Refactor `update_x!` AdditiveNoise

**Files:**
- Modify: `src/largedeviations/sgMAM.jl`

- [ ] **Step 1: Baseline**

```bash
julia --project=. -e 'using Pkg; Pkg.test(; test_args=["sgMAM"])'
```
Expected: PASS.

- [ ] **Step 2: Refactor signature and dispatch**

Replace `update_x!(x, λ, p′, x′′, Hx, ϵ)` with a dispatcher plus an AdditiveNoise method:

```julia
function update_x!(x, λ, p′, x′′, Hx, sys::FreidlinWentzellHamiltonian, ϵ)
    _update_x!(x, λ, p′, x′′, Hx, sys, ϵ)
    return nothing
end

function _update_x!(
        x, λ, p′, x′′, Hx,
        sys::FreidlinWentzellHamiltonian{IIP, D, Hx_t, Hp_t, AF, AdditiveNoise}, ϵ,
    ) where {IIP, D, Hx_t, Hp_t, AF}
    a_diag = LinearAlgebra.diag(sys.a(view(x, :, 1)))
    a_inv  = inv.(a_diag)

    Nx, Nt = size(x); xa = x[:, 1]; xb = x[:, end]; idxc = 2:(Nt - 1)

    Threads.@threads for dof in 1:Nx
        c = a_inv[dof]
        d  = 1 .+ 2 .* ϵ .* λ[2:(end - 1)] .^ 2 .* c
        du = -ϵ .* λ[2:(end - 2)] .^ 2 .* c
        dl = -ϵ .* λ[3:(end - 1)] .^ 2 .* c
        T = LinearAlgebra.Tridiagonal(dl, d, du)
        rhs = Vector{eltype(x)}(undef, length(d))
        @inbounds for k in 1:length(rhs)
            i = k + 1
            rhs[k] = x[dof, i] + ϵ * (λ[i] * p′[dof, i] + Hx[dof, i] - c * λ[i]^2 * x′′[dof, i])
        end
        rhs[1]   += ϵ * c * λ[2]^2 * xa[dof]
        rhs[end] += ϵ * c * λ[end - 1]^2 * xb[dof]
        prob = LinearProblem(T, rhs)
        cache = init(
            prob, LUFactorization();
            alias = SciMLBase.LinearAliasSpecifier(; alias_A = true, alias_b = true),
        )
        solve!(cache)
        x[dof, idxc] .= cache.u
    end
    return nothing
end
```

Note on perf vs the prior main-branch path: the previous implementation built a single `Tridiagonal` once and shared it across the threaded DOF loop. The patch above builds a new `Tridiagonal` per DOF (necessary because anisotropic constant `a` gives different per-DOF coefficients). For the historically-common isotropic case `a ∝ I` all DOFs build the same `Tridiagonal`, which is wasteful. Task 16 introduces per-thread `LinearSolve` caches reused via the cache-update API, which bounds the allocation cost; if the benchmark in Task 25 still shows a regression on the additive isotropic baseline, add a fast-path branch that detects `all(==(a_diag[1]), a_diag)` and reverts to a single shared `Tridiagonal`.

- [ ] **Step 3: Run tests**

```bash
julia --project=. -e 'using Pkg; Pkg.test(; test_args=["sgMAM"])'
```
Expected: PASS (existing tolerances absorb per-DOF reordering noise).

- [ ] **Step 4: Ask for approval, then commit**

```bash
git add src/largedeviations/sgMAM.jl
git commit -m "refactor update_x! to NS-dispatched per-DOF tridiagonal (additive)"
```

---

### Task 13: `update_x!` for `DiagonalNoise`

**Files:**
- Modify: `src/largedeviations/sgMAM.jl`
- Modify: `test/largedeviations/multiplicative_noise.jl`

- [ ] **Step 1: Add an end-to-end test**

Append to `test/largedeviations/multiplicative_noise.jl`:

```julia
@testset "sgMAM end-to-end: 1D OU multiplicative converges" begin
    α = 0.3
    b1d(u, p, t) = SA[-u[1]]
    g1d(u, p, t) = SA[sqrt(1 + α * u[1]^2);;]
    ds = CoupledSDEs(b1d, SA[1.0]; g = g1d, noise_prototype = SMatrix{1, 1}(0.0))
    sys = FreidlinWentzellHamiltonian(ds)

    Nt = 80
    x_initial = reshape(collect(range(1.0, -1.0; length = Nt)), 1, Nt)
    res = minimize_simple_geometric_action(
        sys, x_initial, GeometricGradient(; stepsize = 1.0);
        maxiters = 500, show_progress = false,
    )
    @test isfinite(res.action)
end
```

- [ ] **Step 2: Run, expect failure**

```bash
julia --project=. -e 'using Pkg; Pkg.test(; test_args=["multiplicative_noise"])'
```
Expected: `MethodError` on `DiagonalNoise` `_update_x!`.

- [ ] **Step 3: Implement**

Append to `src/largedeviations/sgMAM.jl`:

```julia
function _update_x!(
        x, λ, p′, x′′, Hx,
        sys::FreidlinWentzellHamiltonian{IIP, D, Hx_t, Hp_t, AF, DiagonalNoise}, ϵ,
    ) where {IIP, D, Hx_t, Hp_t, AF}
    Nx, Nt = size(x)
    a_at = Matrix{eltype(x)}(undef, Nx, Nt)
    for t in 1:Nt
        a_at[:, t] .= LinearAlgebra.diag(sys.a(view(x, :, t)))
    end

    xa = x[:, 1]; xb = x[:, end]; idxc = 2:(Nt - 1)

    Threads.@threads for dof in 1:Nx
        α = [ϵ * λ[i]^2 / a_at[dof, i] for i in 2:(Nt - 1)]
        d  = 1 .+ 2 .* α
        du = -α[1:(end - 1)]
        dl = -α[2:end]
        T = LinearAlgebra.Tridiagonal(dl, d, du)
        rhs = Vector{eltype(x)}(undef, length(d))
        @inbounds for k in 1:length(rhs)
            i = k + 1
            rhs[k] = x[dof, i] + ϵ * (
                λ[i] * p′[dof, i] + Hx[dof, i] - (λ[i]^2 / a_at[dof, i]) * x′′[dof, i]
            )
        end
        rhs[1]   += α[1]   * xa[dof]
        rhs[end] += α[end] * xb[dof]
        prob = LinearProblem(T, rhs)
        cache = init(
            prob, LUFactorization();
            alias = SciMLBase.LinearAliasSpecifier(; alias_A = true, alias_b = true),
        )
        solve!(cache)
        x[dof, idxc] .= cache.u
    end
    return nothing
end
```

- [ ] **Step 4: Run, expect pass**

```bash
julia --project=. -e 'using Pkg; Pkg.test(; test_args=["multiplicative_noise"])'
```
Expected: PASS.

- [ ] **Step 5: Ask for approval, then commit**

```bash
git add src/largedeviations/sgMAM.jl test/largedeviations/multiplicative_noise.jl
git commit -m "add DiagonalNoise update_x! (per-DOF tridiagonal with a(x))"
```

---

### Task 14: `update_x!` for `GeneralNoise` (sparse block-tridiagonal)

**Files:**
- Modify: `src/largedeviations/sgMAM.jl`
- Modify: `src/CriticalTransitions.jl` (add `SparseArrays` import if missing)
- Modify: `test/largedeviations/multiplicative_noise.jl`

- [ ] **Step 1: Add `SparseArrays` import if needed**

```bash
git grep -n 'SparseArrays' -- 'src/'
```
If missing, add `using SparseArrays: SparseArrays, sparse, SparseMatrixCSC` to `src/CriticalTransitions.jl`.

- [ ] **Step 2: Add failing test**

Append to `test/largedeviations/multiplicative_noise.jl`:

```julia
@testset "sgMAM end-to-end: 2D off-diagonal multiplicative converges" begin
    Random.seed!(0)
    b2(u, p, t) = SA[-u[1], -u[2]]
    function g2(u, p, t)
        s11 = 1 + 0.2 * u[1]; s22 = 1 + 0.2 * u[2]
        s12 = 0.3 * u[2];     s21 = 0.3 * u[1]
        return @SMatrix [s11 s12; s21 s22]
    end
    ds = CoupledSDEs(
        b2, SA[1.0, 0.0]; g = g2, noise_prototype = SMatrix{2, 2}(zeros(2, 2)),
    )
    sys = FreidlinWentzellHamiltonian(ds)

    Nt = 60
    xx = collect(range(1.0, 0.0; length = Nt))
    yy = collect(range(0.0, 1.0; length = Nt))
    x_initial = Matrix([xx yy]')
    res = minimize_simple_geometric_action(
        sys, x_initial, GeometricGradient(; stepsize = 1.0);
        maxiters = 300, show_progress = false,
    )
    @test isfinite(res.action)
end
```

- [ ] **Step 3: Run, expect failure**

```bash
julia --project=. -e 'using Pkg; Pkg.test(; test_args=["multiplicative_noise"])'
```
Expected: `MethodError` on `GeneralNoise` `_update_x!`.

- [ ] **Step 4: Implement (uncached first cut)**

Cross-reference: PR #320's `_update_x_general!` (in `gh pr diff 320 --repo JuliaDynamics/CriticalTransitions.jl`, search for the function name) walks the same COO enumeration order. If the indexing math below ever produces a wrong shape, compare against that reference before diving into the math.


Append to `src/largedeviations/sgMAM.jl`:

Derivation (multiply through by `a(x_i)` to avoid matrix inverses):

```
   x_new[i] - ϵ λ_i² a^{-1}(x_i)·(x_new[i-1] - 2 x_new[i] + x_new[i+1]) = R_i
=> (a(x_i) + 2 ϵ λ_i² I) x_new[i] - ϵ λ_i² x_new[i-1] - ϵ λ_i² x_new[i+1]
   = a(x_i) (x[i] + ϵ(λ_i p'_i + Hx_i)) - ϵ λ_i² x''_i
```

So row `i` of the block-tridiagonal carries:
* Diagonal block: `a(x_i) + 2 ϵ λ_i² I`.
* Both off-diagonals from row `i` (to row `i-1` and row `i+1`): `-ϵ λ_i² I`. The crucial subtlety is that block `(i+1, i)` is *row i+1's* left off-diagonal, so it uses `λ_{i+1}²`, **not** `λ_i²`. Enumerate per-row to get this right.

```julia
function _update_x!(
        x, λ, p′, x′′, Hx,
        sys::FreidlinWentzellHamiltonian{IIP, D, Hx_t, Hp_t, AF, GeneralNoise}, ϵ,
    ) where {IIP, D, Hx_t, Hp_t, AF}
    Nx, Nt = size(x)
    N_in   = Nt - 2
    n      = N_in * Nx

    # Precompute a(x_i) at each interior grid point.
    A_blocks = [Matrix(sys.a(view(x, :, i_in + 1))) for i_in in 1:N_in]

    Iv = Int[]; Jv = Int[]; Vv = Vector{eltype(x)}()
    nnz_est = N_in * Nx * Nx + 2 * (N_in - 1) * Nx
    sizehint!(Iv, nnz_est); sizehint!(Jv, nnz_est); sizehint!(Vv, nnz_est)

    rhs = zeros(eltype(x), n)
    xa = x[:, 1]; xb = x[:, end]

    for i_in in 1:N_in
        i = i_in + 1                # full-path index
        rb  = (i_in - 1) * Nx
        A_i = A_blocks[i_in]
        λi2 = λ[i]^2

        # Diagonal block: a(x_i) + 2 ϵ λ_i² I.
        for k1 in 1:Nx, k2 in 1:Nx
            v = A_i[k1, k2] + (k1 == k2 ? 2 * ϵ * λi2 : zero(eltype(x)))
            push!(Iv, rb + k1); push!(Jv, rb + k2); push!(Vv, v)
        end

        # Left off-diagonal block (row i, col i-1): -ϵ λ_i² I.
        if i_in > 1
            base_L = (i_in - 2) * Nx
            for k in 1:Nx
                push!(Iv, rb + k); push!(Jv, base_L + k); push!(Vv, -ϵ * λi2)
            end
        end

        # Right off-diagonal block (row i, col i+1): -ϵ λ_i² I.
        if i_in < N_in
            base_R = i_in * Nx
            for k in 1:Nx
                push!(Iv, rb + k); push!(Jv, base_R + k); push!(Vv, -ϵ * λi2)
            end
        end

        # RHS: A_i (x[i] + ϵ(λ_i p'_i + Hx_i)) - ϵ λ_i² x''_i + boundary contributions.
        rhs_block = A_i * (x[:, i] .+ ϵ .* (λ[i] .* p′[:, i] .+ Hx[:, i])) .-
            ϵ * λi2 .* x′′[:, i]
        if i_in == 1
            rhs_block .+= ϵ * λi2 .* xa
        end
        if i_in == N_in
            rhs_block .+= ϵ * λi2 .* xb
        end
        @inbounds for k in 1:Nx
            rhs[rb + k] = rhs_block[k]
        end
    end

    M = SparseArrays.sparse(Iv, Jv, Vv, n, n)
    prob = LinearProblem(M, rhs)
    cache = init(prob, LUFactorization())
    solve!(cache)
    sol = cache.u
    for i_in in 1:N_in
        rb = (i_in - 1) * Nx
        @inbounds for k in 1:Nx
            x[k, i_in + 1] = sol[rb + k]
        end
    end
    return nothing
end
```

- [ ] **Step 5: Run, expect pass**

```bash
julia --project=. -e 'using Pkg; Pkg.test(; test_args=["multiplicative_noise"])'
```
Expected: PASS.

- [ ] **Step 6: Ask for approval, then commit**

```bash
git add src/largedeviations/sgMAM.jl src/CriticalTransitions.jl test/largedeviations/multiplicative_noise.jl
git commit -m "add GeneralNoise update_x! (sparse block-tridiagonal)"
```

---

### Task 15: Cache sparse pattern and LinearSolve handle (GeneralNoise perf)

**Files:**
- Modify: `src/largedeviations/sgMAM.jl`
- Modify: `test/largedeviations/multiplicative_noise.jl`

- [ ] **Step 1: Add allocation regression test**

Append to `test/largedeviations/multiplicative_noise.jl`:

```julia
@testset "GeneralNoise update_x!: cached pattern allocation regression" begin
    Random.seed!(0)
    b2(u, p, t) = SA[-u[1], -u[2]]
    function g2(u, p, t)
        s11 = 1 + 0.2 * u[1]; s22 = 1 + 0.2 * u[2]
        s12 = 0.3 * u[2];     s21 = 0.3 * u[1]
        return @SMatrix [s11 s12; s21 s22]
    end
    ds = CoupledSDEs(
        b2, SA[1.0, 0.0]; g = g2, noise_prototype = SMatrix{2, 2}(zeros(2, 2)),
    )
    sys = FreidlinWentzellHamiltonian(ds)

    Nt = 40
    xx = collect(range(1.0, 0.0; length = Nt))
    yy = collect(range(0.0, 1.0; length = Nt))
    x = Matrix([xx yy]')
    minimize_simple_geometric_action(
        sys, x, GeometricGradient(; stepsize = 1.0);
        maxiters = 1, show_progress = false,
    )
    bytes_before = Base.gc_num().total_allocd
    minimize_simple_geometric_action(
        sys, x, GeometricGradient(; stepsize = 1.0);
        maxiters = 5, show_progress = false,
    )
    bytes_after = Base.gc_num().total_allocd
    @test (bytes_after - bytes_before) < 5_000_000
end
```

- [ ] **Step 2: Restructure to use a cached workspace**

Build a per-call cache at the top of `minimize_simple_geometric_action(sys::FreidlinWentzellHamiltonian, ...)`:

```julia
struct _SgMAMGeneralCache{T, C}
    M::SparseArrays.SparseMatrixCSC{T, Int}
    idx_map::Vector{Tuple{Int, Int, Int}}  # (block_id, k1, k2) -> nzval index
    off_idx::Vector{Tuple{Int, Int}}       # off-diagonal indices into nzval
    rhs::Vector{T}
    linear_cache::C                        # LinearSolve.LinearCache
end
```

Build the COO arrays once at warmup; convert to `SparseMatrixCSC` via `sparse(I, J, V, n, n)`. Walk the resulting `M.colptr` / `M.rowval` to construct `idx_map`: for each (block, k1, k2) record the index into `M.nzval` where that entry lives. On subsequent iterations, overwrite `M.nzval[idx]` in place for each (block, k1, k2) without rebuilding. Reuse the `LinearSolve.LinearCache` (built once via `init(LinearProblem(M, rhs), LUFactorization())`) and update it via the LinearSolve cache-update API noted in Task 16.

**Threading the cache through the algorithm:**

To avoid breaking the signatures established in Tasks 9 through 13, add the cache as a keyword argument with default `nothing`:

```julia
function update!(x, xdot, xdotdot, p, pdot, lambda, sys, ϵ; cache = nothing)
    central_diff!(pdot, p)
    Hx = sys.H_x(x, p)
    central_diff!(xdotdot, xdot)
    return update_x!(x, lambda, pdot, xdotdot, Hx, sys, ϵ; cache)
end

function update_x!(x, λ, p′, x′′, Hx, sys::FreidlinWentzellHamiltonian, ϵ; cache = nothing)
    _update_x!(x, λ, p′, x′′, Hx, sys, ϵ, cache)
    return nothing
end

# Existing Additive/Diagonal methods ignore the cache argument:
function _update_x!(x, λ, p′, x′′, Hx, sys::...AdditiveNoise, ϵ, _cache)
    # ... existing body ...
end

# General method consumes it:
function _update_x!(x, λ, p′, x′′, Hx, sys::...GeneralNoise, ϵ, cache::_SgMAMGeneralCache)
    # ... cached body using cache.M, cache.idx_map, cache.linear_cache, cache.rhs ...
end
```

In `minimize_geometric_action(sys::FreidlinWentzellHamiltonian, ...)`:

```julia
cache = sys isa FreidlinWentzellHamiltonian{<:Any,<:Any,<:Any,<:Any,<:Any,GeneralNoise} ?
    _build_general_cache(sys, x_initial) :
    nothing
# pass `cache` as the keyword argument to every update!/_sgmam_refresh! call
```

The `nothing` fast path for Additive/Diagonal is type-stable because the method dispatch on `_update_x!` selects the body that ignores `cache` at compile time.

- [ ] **Step 3: Run tests**

```bash
julia --project=. -e 'using Pkg; Pkg.test(; test_args=["multiplicative_noise"])'
```
Expected: PASS.

- [ ] **Step 4: Ask for approval, then commit**

```bash
git add src/largedeviations/sgMAM.jl test/largedeviations/multiplicative_noise.jl
git commit -m "cache sparse pattern and LinearSolve handle in GeneralNoise update_x!"
```

---

### Task 16: Reuse LinearSolve caches in Additive and Diagonal `update_x!` via `reinit!`

**Files:**
- Modify: `src/largedeviations/sgMAM.jl`
- Modify: `src/CriticalTransitions.jl` (add `reinit!` import)

- [ ] **Step 1: Add the import**

In `src/CriticalTransitions.jl`, change:

```julia
using LinearSolve: LinearProblem, LUFactorization, init, solve, solve!
```

to:

```julia
using LinearSolve: LinearSolve, LinearProblem, LUFactorization, init, solve, solve!
```

The cache-update API is `LinearSolve.reinit!(cache; A, b)` in current LinearSolve releases. If the function is not exported under that name in the pinned version, the alternatives in order of preference are:

1. `LinearSolve.set_A!(cache, A)` followed by `LinearSolve.set_b!(cache, b)`.
2. Mutating `cache.A` and `cache.b` directly and calling `LinearSolve.dont_reuse_lu!(cache)` (or the equivalent re-factorization signal) before `solve!`.
3. Rebuilding the `LinearProblem` per call (drops to the Task 14 first-cut perf level).

Pick the working option, verify the test in Step 3 still passes, then update this step's description with the chosen API in the same commit.

- [ ] **Step 2: Refactor Additive and Diagonal `_update_x!` to reuse caches**

For each branch, build `nthreads` caches once at the top of the call with placeholder `Tridiagonal` and `rhs`; in the threaded DOF loop, overwrite the arrays in place and call `LinearSolve.reinit!(cache; A = T, b = rhs)` before `solve!`. Concrete patch for the Additive branch:

```julia
function _update_x!(
        x, λ, p′, x′′, Hx,
        sys::FreidlinWentzellHamiltonian{IIP, D, Hx_t, Hp_t, AF, AdditiveNoise}, ϵ,
    ) where {IIP, D, Hx_t, Hp_t, AF}
    a_diag = LinearAlgebra.diag(sys.a(view(x, :, 1)))
    a_inv  = inv.(a_diag)

    Nx, Nt = size(x); xa = x[:, 1]; xb = x[:, end]; idxc = 2:(Nt - 1)

    nthreads = Threads.nthreads()
    if hasmethod(Threads.nthreads, Tuple{Symbol})
        nthreads = Threads.nthreads(:default) + Threads.nthreads(:interactive)
    end

    caches = map(1:nthreads) do _
        d0  = ones(eltype(x), Nt - 2)
        du0 = zeros(eltype(x), Nt - 3)
        dl0 = zeros(eltype(x), Nt - 3)
        T0  = LinearAlgebra.Tridiagonal(dl0, d0, du0)
        rhs0 = zeros(eltype(x), Nt - 2)
        init(
            LinearProblem(T0, rhs0), LUFactorization();
            alias = SciMLBase.LinearAliasSpecifier(; alias_A = true, alias_b = true),
        )
    end

    Threads.@threads for dof in 1:Nx
        c  = a_inv[dof]
        cache = caches[Threads.threadid()]
        T = cache.A
        @inbounds for i in 1:(Nt - 2)
            T.d[i] = 1 + 2 * ϵ * λ[i + 1]^2 * c
        end
        @inbounds for i in 1:(Nt - 3)
            T.du[i] = -ϵ * λ[i + 1]^2 * c
            T.dl[i] = -ϵ * λ[i + 2]^2 * c
        end
        rhs = cache.b
        @inbounds for k in 1:length(rhs)
            i = k + 1
            rhs[k] = x[dof, i] + ϵ * (λ[i] * p′[dof, i] + Hx[dof, i] - c * λ[i]^2 * x′′[dof, i])
        end
        rhs[1]   += ϵ * c * λ[2]^2 * xa[dof]
        rhs[end] += ϵ * c * λ[end - 1]^2 * xb[dof]
        LinearSolve.reinit!(cache; A = T, b = rhs)
        solve!(cache)
        x[dof, idxc] .= cache.u
    end
    return nothing
end
```

Apply the analogous transform to the `DiagonalNoise` branch.

- [ ] **Step 3: Run tests**

```bash
julia --project=. -e 'using Pkg; Pkg.test(; test_args=["sgMAM", "multiplicative_noise"])'
```
Expected: PASS.

- [ ] **Step 4: Ask for approval, then commit**

```bash
git add src/largedeviations/sgMAM.jl src/CriticalTransitions.jl
git commit -m "reuse LinearSolve caches via reinit! in update_x! (additive, diagonal)"
```

---

## Phase 5: gMAM `geometric_gradient_step!` per `NS`

### Task 17: Add `a_func` and `NS` to `GeometricGradientWorkspace`

**Files:**
- Modify: `src/largedeviations/minimize_geometric_action.jl`
- Modify: `test/largedeviations/gMAM.jl`

- [ ] **Step 1: Add failing test**

Append to `test/largedeviations/gMAM.jl`:

```julia
@testset "GeometricGradientWorkspace carries NoiseShape" begin
    f_lin(u, p, t) = SA[-u[1], -u[2]]
    ds_iso = CoupledSDEs(f_lin, SA[0.0, 0.0]; noise_strength = 1.0)
    xx = collect(range(-1.0, 1.0; length = 20))
    yy = 0.3 .* (-xx .^ 2 .+ 1)
    path = Matrix([xx yy]')
    ws = CT.geometric_gradient_workspace(ds_iso, path)
    @test ws isa CT.GeometricGradientWorkspace{<:Any,<:Any,<:Any,<:Any,<:Any,<:Any,<:Any,<:Any,<:Any,AdditiveNoise}
end
```

- [ ] **Step 2: Run, expect failure**

```bash
julia --project=. -e 'using Pkg; Pkg.test(; test_args=["gMAM"])'
```
Expected: type-parameter mismatch.

- [ ] **Step 3: Extend the workspace**

```julia
struct GeometricGradientWorkspace{Tupdate, Tprime, Tvec, Tmat, Ttmp, Tarc, TA, Tjac, AF, NS <: NoiseShape}
    update::Tupdate
    x_prime::Tprime
    lambdas::Tvec
    lambdas_prime::Tvec
    prod1::Tmat
    prod2::Tmat
    drift_cache::Tmat
    theta_cache::Tmat
    rhs_explicit::Tmat
    temp1::Ttmp
    temp2::Ttmp
    alpha_cache::Tvec
    dl::Tvec
    d::Tvec
    du::Tvec
    v::Tvec
    arc::Tarc
    A::TA
    jac::Tjac
    a_func::AF
end
```

Update `geometric_gradient_workspace` to call `_classify_noise_shape(sys)`, build the trace-normalized `a_func` callable (same logic as in the FreidlinWentzellHamiltonian constructor), and parameterize the workspace by `NS`. For `AdditiveNoise` set `A = inv(a_func(zeros))` (a constant matrix); for state-dependent leave `A = nothing` (consumers always use `ws.a_func`).

- [ ] **Step 4: Run tests**

```bash
julia --project=. -e 'using Pkg; Pkg.test(; test_args=["gMAM"])'
```
Expected: PASS.

- [ ] **Step 5: Ask for approval, then commit**

```bash
git add src/largedeviations/minimize_geometric_action.jl test/largedeviations/gMAM.jl
git commit -m "carry NoiseShape and a_func on GeometricGradientWorkspace"
```

---

### Task 18: Per-`NS` `geometric_gradient_step!` (gMAM diagonal and general)

**Files:**
- Modify: `src/largedeviations/minimize_geometric_action.jl`
- Modify: `test/largedeviations/multiplicative_noise.jl`

- [ ] **Step 1: Add failing tests**

Append to `test/largedeviations/multiplicative_noise.jl`:

```julia
@testset "gMAM diagonal multiplicative converges (vs Adam)" begin
    using OptimizationOptimisers
    Random.seed!(0)
    α = 0.3
    b1d(u, p, t) = SA[-u[1]]
    g1d(u, p, t) = SA[sqrt(1 + α * u[1]^2);;]
    ds = CoupledSDEs(b1d, SA[1.0]; g = g1d, noise_prototype = SMatrix{1, 1}(0.0))

    Nt = 80
    x_initial = reshape(collect(range(1.0, -1.0; length = Nt)), 1, Nt)
    res_g = minimize_geometric_action(
        ds, x_initial, GeometricGradient(; stepsize = 1.0);
        maxiters = 500, show_progress = false,
    )
    res_ref = minimize_geometric_action(
        ds, x_initial, OptimizationOptimisers.Adam(0.01);
        maxiters = 3000, show_progress = false,
    )
    @test isapprox(res_g.action, res_ref.action; rtol = 0.05)
end

@testset "gMAM general multiplicative converges (vs Adam)" begin
    using OptimizationOptimisers
    Random.seed!(0)
    b2(u, p, t) = SA[-u[1], -u[2]]
    function g2(u, p, t)
        s11 = 1 + 0.2 * u[1]; s22 = 1 + 0.2 * u[2]
        s12 = 0.3 * u[2];     s21 = 0.3 * u[1]
        return @SMatrix [s11 s12; s21 s22]
    end
    ds = CoupledSDEs(b2, SA[1.0, 0.0]; g = g2, noise_prototype = SMatrix{2, 2}(zeros(2, 2)))
    Nt = 60
    xx = collect(range(1.0, 0.0; length = Nt))
    yy = collect(range(0.0, 1.0; length = Nt))
    x_initial = Matrix([xx yy]')
    res_g = minimize_geometric_action(
        ds, x_initial, GeometricGradient(; stepsize = 1.0);
        maxiters = 300, show_progress = false,
    )
    res_ref = minimize_geometric_action(
        ds, x_initial, OptimizationOptimisers.Adam(0.01);
        maxiters = 3000, show_progress = false,
    )
    @test isapprox(res_g.action, res_ref.action; rtol = 0.1)
end
```

- [ ] **Step 2: Run, expect failure**

```bash
julia --project=. -e 'using Pkg; Pkg.test(; test_args=["multiplicative_noise"])'
```
Expected: existing `geometric_gradient_step!` calls `covariance_matrix(ds)` which is `nothing` for state-dependent; errors.

- [ ] **Step 3: Implement**

Split the current body into a dispatcher and per-`NS` methods. The implicit step (the shared `Tridiagonal(I - ε λ² ∂_s²)`) is independent of `a` for all three branches per Heymann-Vanden-Eijnden 2008 eq 3.6; only the explicit RHS depends on the noise shape.

```julia
function geometric_gradient_step!(
        ws::GeometricGradientWorkspace{<:Any,<:Any,<:Any,<:Any,<:Any,<:Any,<:Any,<:Any,<:Any,NS},
        sys, path; stepsize = 0.1, diff_order = 4,
    ) where {NS <: NoiseShape}
    _geometric_gradient_step!(ws, sys, path, NS(); stepsize, diff_order)
    return ws.update
end

# Additive branch: existing main-branch logic, reading `A = ws.A` (a cached constant).
function _geometric_gradient_step!(
        ws, sys, path, ::AdditiveNoise; stepsize = 0.1, diff_order = 4,
    )
    N = size(path, 2)
    dx = 1.0 / (N - 1)

    path_velocity!(ws.x_prime, path, ws.arc; order = diff_order)

    @views for i in 2:(N - 1)
        velocity_norm = anorm(ws.x_prime[:, i], ws.A)
        if velocity_norm > 1.0e-14
            ws.drift_cache[:, i - 1] .= drift(sys, path[:, i])
            ws.lambdas[i] = anorm(ws.drift_cache[:, i - 1], ws.A) / velocity_norm
            isfinite(ws.lambdas[i]) || (ws.lambdas[i] = 0.0)
        else
            ws.lambdas[i] = 0.0
        end
    end
    @views for i in 2:(N - 1)
        ws.lambdas_prime[i] = (ws.lambdas[i + 1] - ws.lambdas[i - 1]) / (2 * dx)
    end
    @views for i in 2:(N - 1)
        J = ws.jac(path[:, i])
        LinearAlgebra.mul!(ws.temp1, J, ws.x_prime[:, i])
        LinearAlgebra.mul!(ws.temp2, J', ws.x_prime[:, i])
        ws.prod1[:, i - 1] .= ws.temp1 .- ws.temp2
        LinearAlgebra.mul!(ws.prod2[:, i - 1], J', ws.drift_cache[:, i - 1])
    end

    _gmam_implicit_shared!(ws, path, N, stepsize, dx; is_state_dep = false)
    return nothing
end

# Diagonal multiplicative: per-grid-point a(x_i) (diagonal). λ_i in the a-metric,
# θ_i recovered via diag(a(x_i)), explicit RHS includes (1/2) ∂_l a_kk θ_k² and
# (∂_l a_kk θ_k) φ'_l contributions.
function _geometric_gradient_step!(
        ws, sys, path, ::DiagonalNoise; stepsize = 0.1, diff_order = 4,
    )
    N  = size(path, 2)
    Nx = size(path, 1)
    dx = 1.0 / (N - 1)
    h_fd = max(sqrt(eps(eltype(path))), eltype(path)(1.0e-8))

    path_velocity!(ws.x_prime, path, ws.arc; order = diff_order)

    @views for i in 2:(N - 1)
        ws.drift_cache[:, i - 1] .= drift(sys, path[:, i])
    end

    # λ_i and θ_i (state-dependent diagonal a)
    @views for i in 2:(N - 1)
        xi = path[:, i]
        a_i = LinearAlgebra.diag(ws.a_func(xi))   # diagonal as a vector
        b_i = ws.drift_cache[:, i - 1]
        φp  = ws.x_prime[:, i]
        num = sum(b_i .^ 2 ./ a_i)
        den = sum(φp .^ 2 ./ a_i)
        λ = den > 1.0e-28 ? sqrt(num / den) : 0.0
        isfinite(λ) || (λ = 0.0)
        ws.lambdas[i] = λ
        @inbounds for k in 1:Nx
            ws.theta_cache[k, i - 1] = (λ * φp[k] - b_i[k]) / a_i[k]
        end
    end
    @views for i in 2:(N - 1)
        ws.lambdas_prime[i] = (ws.lambdas[i + 1] - ws.lambdas[i - 1]) / (2 * dx)
    end

    # Explicit RHS at each interior point: -λ H_θφ φ' + a H_φ + λ λ' φ'.
    # H_φ_l = (J^T θ)_l + (1/2) Σ_k θ_k² ∂_l a_kk
    # (H_θφ φ')_k = (J φ')_k + (∂_l a_kk θ_k) φ'_l  (l indexed by k since a is diagonal)
    e = zeros(eltype(path), Nx)
    Hphi  = zeros(eltype(path), Nx)
    Hpphi = zeros(eltype(path), Nx)
    @views for i in 2:(N - 1)
        xi = path[:, i]
        φp = ws.x_prime[:, i]
        θ  = ws.theta_cache[:, i - 1]
        J  = ws.jac(xi)
        fill!(Hphi, 0); fill!(Hpphi, 0)
        LinearAlgebra.mul!(Hpphi, J, φp)
        LinearAlgebra.mul!(Hphi,  J', θ)
        @inbounds for l in 1:Nx
            fill!(e, 0); e[l] = h_fd
            ap = LinearAlgebra.diag(ws.a_func(xi .+ e))
            am = LinearAlgebra.diag(ws.a_func(xi .- e))
            @inbounds for k in 1:Nx
                dla = (ap[k] - am[k]) / (2 * h_fd)
                Hphi[l] += 0.5 * dla * θ[k]^2
                Hpphi[k] += dla * θ[k] * φp[l]
            end
        end
        a_i = LinearAlgebra.diag(ws.a_func(xi))
        λ   = ws.lambdas[i]
        λp  = ws.lambdas_prime[i]
        @inbounds for k in 1:Nx
            ws.rhs_explicit[k, i - 1] = -λ * Hpphi[k] + a_i[k] * Hphi[k] + λ * λp * φp[k]
        end
    end

    _gmam_implicit_shared!(ws, path, N, stepsize, dx; is_state_dep = true)
    return nothing
end

# General multiplicative: per-grid-point matrix a(x_i). θ_i = a(x_i)^{-1} (λ φ' - b).
# Explicit RHS uses (1/2) θ^T (∂_l a) θ for each l, plus (∂_l a · θ) φ'_l terms.
function _geometric_gradient_step!(
        ws, sys, path, ::GeneralNoise; stepsize = 0.1, diff_order = 4,
    )
    N  = size(path, 2)
    Nx = size(path, 1)
    dx = 1.0 / (N - 1)
    h_fd = max(sqrt(eps(eltype(path))), eltype(path)(1.0e-8))

    path_velocity!(ws.x_prime, path, ws.arc; order = diff_order)

    @views for i in 2:(N - 1)
        ws.drift_cache[:, i - 1] .= drift(sys, path[:, i])
    end

    @views for i in 2:(N - 1)
        xi = path[:, i]
        a_i = ws.a_func(xi)
        b_i = ws.drift_cache[:, i - 1]
        φp  = ws.x_prime[:, i]
        A_inv_b   = a_i \ Vector(b_i)
        A_inv_phi = a_i \ Vector(φp)
        num = dot(b_i, A_inv_b)
        den = dot(φp, A_inv_phi)
        λ = den > 1.0e-28 ? sqrt(num / den) : 0.0
        isfinite(λ) || (λ = 0.0)
        ws.lambdas[i] = λ
        θ = a_i \ (λ .* Vector(φp) .- Vector(b_i))
        @inbounds for k in 1:Nx
            ws.theta_cache[k, i - 1] = θ[k]
        end
    end
    @views for i in 2:(N - 1)
        ws.lambdas_prime[i] = (ws.lambdas[i + 1] - ws.lambdas[i - 1]) / (2 * dx)
    end

    e = zeros(eltype(path), Nx)
    Hphi  = zeros(eltype(path), Nx)
    Hpphi = zeros(eltype(path), Nx)
    aHphi = zeros(eltype(path), Nx)
    @views for i in 2:(N - 1)
        xi = path[:, i]
        φp = ws.x_prime[:, i]
        θ  = ws.theta_cache[:, i - 1]
        J  = ws.jac(xi)
        a_i = ws.a_func(xi)
        fill!(Hphi, 0); fill!(Hpphi, 0)
        LinearAlgebra.mul!(Hpphi, J, φp)
        LinearAlgebra.mul!(Hphi,  J', θ)
        @inbounds for l in 1:Nx
            fill!(e, 0); e[l] = h_fd
            ap = ws.a_func(xi .+ e)
            am = ws.a_func(xi .- e)
            @inbounds for k1 in 1:Nx, k2 in 1:Nx
                dla = (ap[k1, k2] - am[k1, k2]) / (2 * h_fd)
                Hphi[l] += 0.5 * dla * θ[k1] * θ[k2]
                Hpphi[k1] += dla * θ[k2] * φp[l]
            end
        end
        LinearAlgebra.mul!(aHphi, a_i, Hphi)
        λ  = ws.lambdas[i]; λp = ws.lambdas_prime[i]
        @inbounds for k in 1:Nx
            ws.rhs_explicit[k, i - 1] = -λ * Hpphi[k] + aHphi[k] + λ * λp * φp[k]
        end
    end

    _gmam_implicit_shared!(ws, path, N, stepsize, dx; is_state_dep = true)
    return nothing
end

# Implicit step: (I - ε λ_i² ∂_s² / dx²) acting per DOF; shared across noise shapes.
# Matches the existing main-branch gMAM implementation: row i carries `α_i = ε λ_i²/dx²`,
# and BOTH off-diagonals from row i use `α_i` (not the neighbors' values). The
# `α_i == 0` skip mirrors the main code's behavior when a path point has λ=0.
function _gmam_implicit_shared!(ws, path, N, stepsize, dx; is_state_dep)
    ws.d  .= 1
    ws.dl .= 0
    ws.du .= 0
    @views for i in 2:(N - 1)
        ws.alpha_cache[i] = stepsize * ws.lambdas[i]^2 / dx^2
        if isfinite(ws.alpha_cache[i])
            ws.d[i]      += 2 * ws.alpha_cache[i]
            ws.dl[i - 1]  = -ws.alpha_cache[i]
            ws.du[i]      = -ws.alpha_cache[i]
        else
            ws.alpha_cache[i] = 0.0
        end
    end

    T = LinearAlgebra.Tridiagonal(ws.dl, ws.d, ws.du)
    prob = LinearProblem(T, ws.v)
    cache = init(
        prob, LUFactorization();
        alias = SciMLBase.LinearAliasSpecifier(; alias_A = true, alias_b = true),
    )
    rhs = cache.b

    @inbounds for j in 1:size(path, 1)
        rhs[1]   = path[j, 1]
        rhs[end] = path[j, end]
        for i in 2:(N - 1)
            if ws.alpha_cache[i] == 0.0
                rhs[i] = path[j, i]
                continue
            end
            rhs_val = if is_state_dep
                path[j, i] + stepsize * ws.rhs_explicit[j, i - 1]
            else
                path[j, i] + stepsize * (
                    -ws.lambdas[i] * ws.prod1[j, i - 1] - ws.prod2[j, i - 1] +
                        ws.lambdas[i] * ws.lambdas_prime[i] * ws.x_prime[j, i]
                )
            end
            rhs[i] = isfinite(rhs_val) ? rhs_val : path[j, i]
        end
        solve!(cache)
        ws.update[j, :] .= cache.u
    end
    return nothing
end
```

The `ws.v` buffer (length `N`) is reused as the LinearProblem RHS storage; the existing `geometric_gradient_workspace` already allocates it.

Notes on adaptation from PR #320:
* The PR #320 version uses a `LinearProblem` per DOF inside the implicit step; the patch above keeps that pattern (rebuild `T` once, share across DOFs since it does not depend on the DOF index). Caching the LinearSolve handle across iterations is the work of a follow-up if benchmarks show it matters.
* PR #320 has a separate `_normalized_a_shape_func` for trace normalization. Here the workspace's `ws.a_func` is already trace-normalized at construction (Task 17), so we evaluate it directly.

- [ ] **Step 4: Run tests**

```bash
julia --project=. -e 'using Pkg; Pkg.test(; test_args=["gMAM", "multiplicative_noise"])'
```
Expected: PASS.

- [ ] **Step 5: Ask for approval, then commit**

```bash
git add src/largedeviations/minimize_geometric_action.jl test/largedeviations/multiplicative_noise.jl
git commit -m "extend geometric_gradient_step! to diagonal and general multiplicative noise"
```

---

## Phase 6: Action functions

### Task 19: Add `_action_metric`

**Files:**
- Modify: `src/largedeviations/action.jl`
- Modify: `test/largedeviations/multiplicative_noise.jl`

- [ ] **Step 1: Add failing test**

Append:

```julia
@testset "_action_metric: additive returns constant inv" begin
    f_lin(u, p, t) = SA[-u[1], -u[2]]
    ds = CoupledSDEs(f_lin, SA[0.0, 0.0]; noise_strength = 2.0)
    metric = CT._action_metric(ds)
    @test metric isa Base.Returns
    @test metric(zeros(2)) ≈ inv(LinearAlgebra.Diagonal([1.0, 1.0]))
end

@testset "_action_metric: state-dependent evaluates per point" begin
    α = 0.3
    b1d(u, p, t) = SA[-u[1]]
    g1d(u, p, t) = SA[sqrt(1 + α * u[1]^2);;]
    ds = CoupledSDEs(b1d, SA[1.0]; g = g1d, noise_prototype = SMatrix{1, 1}(0.0))
    metric = CT._action_metric(ds)
    @test !(metric isa Base.Returns)
    @test metric([1.0])[1, 1] ≈ 1.0
end
```

- [ ] **Step 2: Run, expect failure**

```bash
julia --project=. -e 'using Pkg; Pkg.test(; test_args=["multiplicative_noise"])'
```
Expected: `UndefVarError: _action_metric not defined`.

- [ ] **Step 3: Implement**

Add to `src/largedeviations/action.jl` (top, after the imports):

```julia
"""
    _action_metric(sys::CoupledSDEs)

Returns `x -> A(x)`, where `A(x)` is the inverse of the trace-normalized diffusion
tensor `a(x)/s` with `s = tr(a(u₀))/D`. For additive `sys` the callable is `Returns(A)`
(constant); for state-dependent it evaluates `σ(x)σ(x)ᵀ / s` per call.
"""
function _action_metric(sys::CoupledSDEs)
    _classify_noise_shape(sys)
    σ_fn = diffusion_function(sys)
    u₀ = current_state(sys); ps = current_parameters(sys)
    a_of(x) = let σx = σ_fn(x, ps, 0.0)
        σ_mat = σx isa AbstractMatrix ? σx : LinearAlgebra.Diagonal(σx)
        σ_mat * σ_mat'
    end
    a0 = a_of(u₀)
    s  = LinearAlgebra.tr(a0) / size(a0, 1)
    if sys.noise_type[:additive]
        return Returns(inv(a0 / s))
    else
        return x -> inv(a_of(x) / s)
    end
end
```

- [ ] **Step 4: Run, expect pass**

```bash
julia --project=. -e 'using Pkg; Pkg.test(; test_args=["multiplicative_noise"])'
```
Expected: PASS.

- [ ] **Step 5: Ask for approval, then commit**

```bash
git add src/largedeviations/action.jl test/largedeviations/multiplicative_noise.jl
git commit -m "add _action_metric helper for FW action integrand"
```

---

### Task 20: Generalize `fw_action`, `om_action`, `geometric_action`

**Files:**
- Modify: `src/largedeviations/action.jl`
- Modify: `test/largedeviations/multiplicative_noise.jl`

- [ ] **Step 1: Add failing analytic test**

Append:

```julia
@testset "fw_action: 1D OU multiplicative vs analytic Simpson" begin
    α = 0.3
    b1d(u, p, t) = SA[-u[1]]
    g1d(u, p, t) = SA[sqrt(1 + α * u[1]^2);;]
    ds = CoupledSDEs(b1d, SA[1.0]; g = g1d, noise_prototype = SMatrix{1, 1}(0.0))

    N, T = 200, 1.0
    path = reduce(hcat, range([1.0], [0.0]; length = N))
    time = range(0.0, T; length = N)
    S = fw_action(ds, path, time)

    function simpson(f, a, b, n)
        h = (b - a) / n; s = f(a) + f(b)
        for i in 1:2:(n - 1); s += 4 * f(a + i * h); end
        for i in 2:2:(n - 2); s += 2 * f(a + i * h); end
        return s * h / 3
    end
    s_norm = 1 + α
    integrand = t -> t^2 / (1 + α * (1 - t)^2)
    analytic = s_norm * simpson(integrand, 0, 1, 10_000) / 2
    @test isapprox(S, analytic; rtol = 1e-4)
end
```

- [ ] **Step 2: Run, expect failure**

```bash
julia --project=. -e 'using Pkg; Pkg.test(; test_args=["multiplicative_noise"])'
```
Expected: existing `fw_action` errors because `covariance_matrix(sys) === nothing` for state-dependent.

- [ ] **Step 3: Implement**

In `src/largedeviations/action.jl`, refactor:

```julia
function fw_action(sys::CoupledSDEs, path, time)
    @assert all(diff(time) .≈ diff(time[1:2])) "Freidlin-Wentzell action is only defined for equispaced time"
    A_at = _action_metric(sys)
    integrand = fw_integrand(sys, path, time, A_at)
    S = 0
    for i in 1:(size(path, 2) - 1)
        S += (integrand[i + 1] + integrand[i]) / 2 * (time[i + 1] - time[i])
    end
    return S / 2
end

function fw_integrand(sys::CoupledSDEs, path, time, A_at)
    v = path_velocity(path, time; order = 4)
    sqnorm = zeros(size(path, 2))
    b(x) = drift(sys, x)
    for i in axes(path, 2)
        drift_i = b(path[:, i])
        A_i = A_at(path[:, i])
        sqnorm[i] = anorm(v[:, i] - drift_i, A_i; square = true)
    end
    return sqnorm
end
```

Apply the same `A_at = _action_metric(sys)` plus per-grid-point evaluation to `om_action` and `geometric_action(sys::CoupledSDEs, path, arclength)`. The drift-only overload `geometric_action(b::Function, path, arclength; A)` keeps its existing signature; users wanting a constant metric pass an `AbstractMatrix` for `A`.

- [ ] **Step 4: Run full suite**

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```
Expected: PASS. Existing additive tests use `Returns(A)` and collapse to the constant path.

- [ ] **Step 5: Ask for approval, then commit**

```bash
git add src/largedeviations/action.jl test/largedeviations/multiplicative_noise.jl
git commit -m "generalize fw_action, om_action, geometric_action to state-dependent a"
```

---

## Phase 7: Unified entry point and #263

### Task 21: Add `minimize_geometric_action(::FreidlinWentzellHamiltonian, ...)`

**Files:**
- Modify: `src/largedeviations/minimize_geometric_action.jl`
- Create: `test/largedeviations/unified_api.jl`
- Modify: `test/runtests.jl`

- [ ] **Step 1: Create the test file**

```julia
using CriticalTransitions, StaticArrays
using Test
using LinearAlgebra

const CT = CriticalTransitions

@testset "minimize_geometric_action dispatches: CoupledSDEs (gMAM) ≈ FreidlinWentzellHamiltonian (sgMAM)" begin
    function meier_stein(u, p, t)
        x, y = u
        return SA[x - x^3 - 10 * x * y^2, -(1 + x^2) * y]
    end
    σ = 0.25
    ds = CoupledSDEs(meier_stein, zeros(2); noise_strength = σ)
    sys_h = FreidlinWentzellHamiltonian(ds)

    Nt = 60
    xx = range(-1.0, 1.0; length = Nt)
    yy = 0.3 .* (-xx .^ 2 .+ 1)
    x_initial = Matrix([xx yy]')

    res_g = minimize_geometric_action(
        ds, x_initial, GeometricGradient(; stepsize = 1.0);
        maxiters = 500, show_progress = false,
    )
    res_s = minimize_geometric_action(
        sys_h, x_initial, GeometricGradient(; stepsize = 1.0);
        maxiters = 500, show_progress = false,
    )
    @test isapprox(res_g.action, res_s.action; rtol = 1e-2)
end
```

Add `include("largedeviations/unified_api.jl")` to `test/runtests.jl`.

- [ ] **Step 2: Run, expect failure**

```bash
julia --project=. -e 'using Pkg; Pkg.test(; test_args=["unified_api"])'
```
Expected: `MethodError` on `minimize_geometric_action(::FreidlinWentzellHamiltonian, ...)`.

- [ ] **Step 3: Implement**

In `src/largedeviations/minimize_geometric_action.jl`:

```julia
function minimize_geometric_action(
        sys::FreidlinWentzellHamiltonian,
        x_initial::Matrix,
        optimizer = GeometricGradient(; stepsize = 1.0e3);
        kwargs...,
    )
    return minimize_simple_geometric_action(sys, x_initial, optimizer; kwargs...)
end

function minimize_geometric_action(
        sys::FreidlinWentzellHamiltonian, x_initial::StateSpaceSet, optimizer; kwargs...,
    )
    return minimize_simple_geometric_action(sys, x_initial, optimizer; kwargs...)
end
```

- [ ] **Step 4: Run, expect pass**

```bash
julia --project=. -e 'using Pkg; Pkg.test(; test_args=["unified_api"])'
```
Expected: PASS.

- [ ] **Step 5: Ask for approval, then commit**

```bash
git add src/largedeviations/minimize_geometric_action.jl test/largedeviations/unified_api.jl test/runtests.jl
git commit -m "add minimize_geometric_action dispatch for FreidlinWentzellHamiltonian"
```

---

### Task 22: Remove `minimize_simple_geometric_action` and tighten `::CoupledSDEs` (#263)

**Files:**
- Modify: `src/largedeviations/sgMAM.jl`
- Modify: `src/largedeviations/minimize_geometric_action.jl`
- Modify: `src/CriticalTransitions.jl`
- Modify: `test/largedeviations/sgMAM.jl`
- Modify: `test/largedeviations/unified_api.jl`
- Modify: examples and benchmarks

- [ ] **Step 1: Add failing rejection test**

Append to `test/largedeviations/unified_api.jl`:

```julia
@testset "minimize_geometric_action rejects CoupledODEs (#263)" begin
    f_lin(u, p, t) = SA[-u[1], -u[2]]
    ds_ode = CoupledODEs(f_lin, SA[0.0, 0.0])
    xx = collect(range(-1.0, 1.0; length = 20))
    yy = 0.3 .* (-xx .^ 2 .+ 1)
    path = Matrix([xx yy]')
    @test_throws MethodError minimize_geometric_action(ds_ode, path)
end

@testset "minimize_simple_geometric_action is removed (#326)" begin
    @test !isdefined(CriticalTransitions, :minimize_simple_geometric_action)
end
```

- [ ] **Step 2: Run, expect failure**

```bash
julia --project=. -e 'using Pkg; Pkg.test(; test_args=["unified_api"])'
```
Expected: `MethodError` not thrown for CoupledODEs (existing `::ContinuousTimeDynamicalSystem` dispatch covers it); `isdefined` returns true.

- [ ] **Step 3: Make the changes**

1. In `src/largedeviations/minimize_geometric_action.jl`, change every `sys::ContinuousTimeDynamicalSystem` to `sys::CoupledSDEs` (four `minimize_geometric_action` overloads, `_gmam_setup`, `geometric_gradient_workspace`, `geometric_gradient_step`).

2. Inline `minimize_simple_geometric_action`'s body into `minimize_geometric_action(::FreidlinWentzellHamiltonian, ...)` and delete the `minimize_simple_geometric_action` overloads in `src/largedeviations/sgMAM.jl`.

3. In `src/CriticalTransitions.jl`, change `export minimize_simple_geometric_action, FreidlinWentzellHamiltonian` to `export FreidlinWentzellHamiltonian`.

4. Update **every test file that references `minimize_simple_geometric_action`**: rename to `minimize_geometric_action`. To enumerate:
   ```bash
   git grep -nl 'minimize_simple_geometric_action' -- 'test/'
   ```
   Expect at least `test/largedeviations/sgMAM.jl` and `test/largedeviations/multiplicative_noise.jl` (the latter was written in Tasks 10, 13, 14, 15 using the soon-to-be-removed name).

5. Update `examples/sgMAM_KPO.jl`, `examples/backtracking_KPO.jl`, and any benchmark scripts the same way.

- [ ] **Step 4: Run full suite**

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```
Expected: PASS.

- [ ] **Step 5: Ask for approval, then commit**

```bash
git status
git add -u
git commit -m "remove minimize_simple_geometric_action; tighten gMAM to ::CoupledSDEs"
```

---

## Phase 8: Performance verification

### Task 23: `@inferred` type-stability tests

**Files:**
- Modify: `test/largedeviations/sgMAM.jl`
- Modify: `test/largedeviations/gMAM.jl`

- [ ] **Step 1: Add the tests**

Append to `test/largedeviations/sgMAM.jl`:

```julia
@testset "FreidlinWentzellHamiltonian inner-loop type-stable" begin
    function meier_stein(u, p, t)
        x, y = u
        return SA[x - x^3 - 10 * x * y^2, -(1 + x^2) * y]
    end
    ds = CoupledSDEs(meier_stein, zeros(2); noise_strength = 0.25)
    sys = FreidlinWentzellHamiltonian(ds)

    Nt = 20
    xx = range(-1.0, 1.0; length = Nt)
    yy = 0.3 .* (-xx .^ 2 .+ 1)
    x = Matrix([xx yy]')
    xdot = zeros(size(x)); p = zeros(size(x)); λ = zeros(1, Nt)
    pdot = zeros(size(x)); xdotdot = zeros(size(x))
    CT.central_diff!(xdot, x)
    CT.update_p!(p, λ, x, xdot, sys)
    Hx = sys.H_x(x, p)
    CT.central_diff!(pdot, p); CT.central_diff!(xdotdot, xdot)

    @inferred CT.update_p!(p, λ, x, xdot, sys)
    @inferred CT.update_x!(x, λ, pdot, xdotdot, Hx, sys, 1.0)
end
```

Append to `test/largedeviations/gMAM.jl`:

```julia
@testset "GeometricGradientWorkspace + step! type-stable" begin
    function meier_stein(u, p, t)
        x, y = u
        return SA[x - x^3 - 10 * x * y^2, -(1 + x^2) * y]
    end
    ds = CoupledSDEs(meier_stein, zeros(2); noise_strength = 0.25)
    Nt = 20
    xx = range(-1.0, 1.0; length = Nt)
    yy = 0.3 .* (-xx .^ 2 .+ 1)
    path = Matrix([xx yy]')
    ws = CT.geometric_gradient_workspace(ds, path)
    @inferred CT.geometric_gradient_step!(ws, ds, path; stepsize = 0.1)
end
```

- [ ] **Step 2: Run; if `@inferred` fails, fix instabilities**

```bash
julia --project=. -e 'using Pkg; Pkg.test(; test_args=["sgMAM", "gMAM"])'
```
Expected: PASS. Common culprits: the `a_callable` closure not returning a concretely-typed matrix (ensure `Returns(...)` wraps `Diagonal{Float64, Vector{Float64}}`), or `_action_metric` returning a `Union{Returns, Function}` (split via dispatch instead of a runtime branch).

- [ ] **Step 3: Ask for approval, then commit**

```bash
git add test/largedeviations/sgMAM.jl test/largedeviations/gMAM.jl
git commit -m "type-stability regression tests for hot paths"
```

---

### Task 24: Allocation regression test (Additive)

**Files:**
- Modify: `test/largedeviations/sgMAM.jl`

- [ ] **Step 1: Add the test**

Append:

```julia
@testset "AdditiveNoise update_x! no per-iteration sparse alloc" begin
    function meier_stein(u, p, t)
        x, y = u
        return SA[x - x^3 - 10 * x * y^2, -(1 + x^2) * y]
    end
    ds = CoupledSDEs(meier_stein, zeros(2); noise_strength = 0.25)
    sys = FreidlinWentzellHamiltonian(ds)
    Nt = 60
    xx = range(-1.0, 1.0; length = Nt)
    yy = 0.3 .* (-xx .^ 2 .+ 1)
    x_initial = Matrix([xx yy]')
    minimize_geometric_action(
        sys, x_initial, GeometricGradient(; stepsize = 1.0);
        maxiters = 1, show_progress = false,
    )
    bytes_before = Base.gc_num().total_allocd
    minimize_geometric_action(
        sys, x_initial, GeometricGradient(; stepsize = 1.0);
        maxiters = 10, show_progress = false,
    )
    bytes_after = Base.gc_num().total_allocd
    @test (bytes_after - bytes_before) < 5_000_000
end
```

- [ ] **Step 2: Run**

```bash
julia --project=. -e 'using Pkg; Pkg.test(; test_args=["sgMAM"])'
```
Expected: PASS.

- [ ] **Step 3: Ask for approval, then commit**

```bash
git add test/largedeviations/sgMAM.jl
git commit -m "allocation regression test for additive update_x!"
```

---

### Task 25: Add `benchmarks/multiplicative_noise.jl`

**Files:**
- Create: `benchmarks/multiplicative_noise.jl`

- [ ] **Step 1: Create the benchmark file**

```julia
using CriticalTransitions, StaticArrays, BenchmarkTools
using LinearAlgebra
using Random

const SUITE = BenchmarkGroup()

let
    function meier_stein(u, p, t)
        x, y = u
        return SA[x - x^3 - 10 * x * y^2, -(1 + x^2) * y]
    end
    ds = CoupledSDEs(meier_stein, zeros(2); noise_strength = 0.25)
    Nt = 100
    xx = range(-1.0, 1.0; length = Nt)
    yy = 0.3 .* (-xx .^ 2 .+ 1)
    x_initial = Matrix([xx yy]')
    SUITE["sgMAM/additive/meier_stein"] = @benchmarkable begin
        minimize_geometric_action(
            $(FreidlinWentzellHamiltonian(ds)), $x_initial,
            GeometricGradient(; stepsize = 1.0);
            maxiters = 100, show_progress = false,
        )
    end
end

let
    α = 0.3
    b1d(u, p, t) = SA[-u[1]]
    g1d(u, p, t) = SA[sqrt(1 + α * u[1]^2);;]
    ds = CoupledSDEs(b1d, SA[1.0]; g = g1d, noise_prototype = SMatrix{1, 1}(0.0))
    Nt = 80
    x_initial = reshape(collect(range(1.0, -1.0; length = Nt)), 1, Nt)
    SUITE["sgMAM/diagonal/1d_ou"] = @benchmarkable begin
        minimize_geometric_action(
            $(FreidlinWentzellHamiltonian(ds)), $x_initial,
            GeometricGradient(; stepsize = 1.0);
            maxiters = 100, show_progress = false,
        )
    end
end

let
    b2(u, p, t) = SA[-u[1], -u[2]]
    function g2(u, p, t)
        s11 = 1 + 0.2 * u[1]; s22 = 1 + 0.2 * u[2]
        s12 = 0.3 * u[2];     s21 = 0.3 * u[1]
        return @SMatrix [s11 s12; s21 s22]
    end
    ds = CoupledSDEs(b2, SA[1.0, 0.0]; g = g2, noise_prototype = SMatrix{2, 2}(zeros(2, 2)))
    Nt = 60
    xx = collect(range(1.0, 0.0; length = Nt))
    yy = collect(range(0.0, 1.0; length = Nt))
    x_initial = Matrix([xx yy]')
    SUITE["sgMAM/general/2d_offdiag"] = @benchmarkable begin
        minimize_geometric_action(
            $(FreidlinWentzellHamiltonian(ds)), $x_initial,
            GeometricGradient(; stepsize = 1.0);
            maxiters = 100, show_progress = false,
        )
    end
end

SUITE
```

- [ ] **Step 2: Sanity-check load**

```bash
julia --project=benchmarks -e 'include("benchmarks/multiplicative_noise.jl"); println(keys(SUITE))'
```
Expected: three benchmark keys print.

- [ ] **Step 3: Run additive baseline against `main` (perf-gate the PR per spec section 3.4)**

The spec commits to the additive-isotropic sgMAM path staying within `1.05×` of `main`'s descent cost. Verify before declaring the implementation done:

```bash
# Run the new branch's Meier-Stein baseline:
julia --project=benchmarks -e '
    include("benchmarks/multiplicative_noise.jl");
    r = run(SUITE["sgMAM/additive/meier_stein"]; samples = 10, seconds = 30);
    println(minimum(r));
'

# Then check out main, run the same harness, and compare:
git stash push
git checkout main
julia --project=benchmarks -e '
    include("benchmarks/kpo.jl");  # or the closest equivalent on main
    # ... record baseline ...
'
git checkout fw-hamiltonian-multiplicative-noise
git stash pop
```

If the new branch's `minimum(r).time` exceeds `1.05 × baseline.time`, the PR holds. Mitigation options (in order): add an isotropic-detection fast path (see Task 12 note), enable a per-thread LinearSolve cache via Task 16's `reinit!`, or accept the regression with a written justification.

- [ ] **Step 4: Ask for approval, then commit**

```bash
git add benchmarks/multiplicative_noise.jl
git commit -m "add benchmarks for multiplicative noise paths"
```

---

## Phase 9: Remaining tests and docs

### Task 26: Additive non-diagonal Σ test

**Files:**
- Modify: `test/largedeviations/multiplicative_noise.jl`

- [ ] **Step 1: Add the test**

```julia
@testset "Additive non-diagonal Σ goes through GeneralNoise" begin
    function meier_stein(u, p, t)
        x, y = u
        return SA[x - x^3 - 10 * x * y^2, -(1 + x^2) * y]
    end
    Nt = 60
    xx = range(-1.0, 1.0; length = Nt)
    yy = 0.3 .* (-xx .^ 2 .+ 1)
    x_initial = Matrix([xx yy]')

    θ = 0.4
    R = [cos(θ) -sin(θ); sin(θ) cos(θ)]
    D = Diagonal([0.5, 2.0])
    Q_rot = R * D * R'
    ds_rot = CoupledSDEs(meier_stein, zeros(2); covariance = Q_rot)
    sys_rot = FreidlinWentzellHamiltonian(ds_rot)
    @test sys_rot isa FreidlinWentzellHamiltonian{<:Any, 2, <:Any, <:Any, <:Any, GeneralNoise}

    res = minimize_geometric_action(
        sys_rot, x_initial, GeometricGradient(; stepsize = 1.0);
        maxiters = 500, show_progress = false,
    )
    @test isfinite(res.action)
end
```

- [ ] **Step 2: Run**

```bash
julia --project=. -e 'using Pkg; Pkg.test(; test_args=["multiplicative_noise"])'
```
Expected: PASS.

- [ ] **Step 3: Ask for approval, then commit**

```bash
git add test/largedeviations/multiplicative_noise.jl
git commit -m "add additive non-diagonal Σ regression"
```

---

### Task 27: Rank-deficient rejection test

**Files:**
- Create: `test/largedeviations/rank_deficient_rejection.jl`
- Modify: `test/runtests.jl`

- [ ] **Step 1: Create the test file**

```julia
using CriticalTransitions, StaticArrays
using Test

const CT = CriticalTransitions

@testset "Rank-deficient rejection (#325 deferred)" begin
    function langevin(u, p, t)
        x, p_ = u
        return SA[p_, -x - 0.1 * p_]
    end
    g_langevin(u, p, t) = SA[0.0 0.0; 0.0 sqrt(0.2)]
    ds = CoupledSDEs(
        langevin, SA[0.0, 0.0]; g = g_langevin,
        noise_prototype = SMatrix{2, 2}(zeros(2, 2)),
    )

    err_h = try
        FreidlinWentzellHamiltonian(ds); nothing
    catch e
        e
    end
    @test err_h isa ArgumentError
    @test occursin("rank-deficient", err_h.msg)
    # Message names a workaround (FreidlinWentzellHamiltonian hand-rolled constructor).
    @test occursin("FreidlinWentzellHamiltonian", err_h.msg)

    Nt = 20
    xx = range(-1.0, 1.0; length = Nt)
    yy = 0.3 .* (-xx .^ 2 .+ 1)
    path = Matrix([xx yy]')
    err_g = try
        minimize_geometric_action(ds, path); nothing
    catch e
        e
    end
    @test err_g isa ArgumentError
    @test occursin("rank-deficient", err_g.msg)
end
```

Add `include("largedeviations/rank_deficient_rejection.jl")` to `test/runtests.jl`.

- [ ] **Step 2: Run, expect pass**

```bash
julia --project=. -e 'using Pkg; Pkg.test(; test_args=["rank_deficient_rejection"])'
```
Expected: PASS.

- [ ] **Step 3: Ask for approval, then commit**

```bash
git add test/largedeviations/rank_deficient_rejection.jl test/runtests.jl
git commit -m "add rank-deficient rejection regression tests"
```

---

### Task 28: Docs: Hamiltonian picture, method table, NoiseShape

**Files:**
- Modify: `docs/src/man/largedeviations.md`
- Modify: `docs/src/examples/sgMAM_KPO.md` (if it exists)

- [ ] **Step 1: Update `docs/src/man/largedeviations.md`**

1. Add a "Hamiltonian picture" subsection after the FW intro:

````markdown
## Hamiltonian picture

The Freidlin-Wentzell rate function admits an equivalent Hamiltonian formulation with
conjugate momentum `p`:

```math
H(\varphi, p) = \langle b(\varphi), p \rangle + \tfrac{1}{2}\, \langle p, a(\varphi)\, p \rangle,
```

where `a(x) = σ(x)σ(x)ᵀ`. The action along an instanton (where `H ≡ 0`) reduces to
`S[φ] = ∫_0^T ⟨p, dφ⟩ dt`, the form used by `minimize_geometric_action` when the input
is a `FreidlinWentzellHamiltonian`.

The implementation classifies `a(x)` once at construction into one of three
`NoiseShape` singletons (`AdditiveNoise`, `DiagonalNoise`, `GeneralNoise`). The tag is a
type parameter, so inner-loop dispatch happens at compile time.

Rank-deficient `a(x)` is not supported in this release; see issue #325.
````

2. Update the method table to reflect: sgMAM now supports state-dependent and general noise; rank-deficient is the only remaining limitation.

3. Replace all `ExtendedPhaseSpace` references with `FreidlinWentzellHamiltonian`.

4. Remove the "current sgMAM implementation assumes additive noise" caveat.

- [ ] **Step 2: Build the docs**

```bash
julia --project=docs -e 'using Documenter, CriticalTransitions; include("docs/make.jl")'
```
Expected: docs build without errors.

- [ ] **Step 3: Ask for approval, then commit**

```bash
git add docs/src/man/largedeviations.md docs/src/examples/sgMAM_KPO.md
git commit -m "document Hamiltonian picture and NoiseShape dispatch"
```

---

### Task 29: Final cleanup

**Files:**
- Various

- [ ] **Step 1: Find leftover references**

```bash
git grep -n 'ExtendedPhaseSpace\|minimize_simple_geometric_action\|proper_sgMAM_system' -- ':!docs/superpowers/'
```
Expected: empty.

- [ ] **Step 2: Check `test/code_quality.jl`**

If it lists `ExtendedPhaseSpace` in any exclusion set, remove that entry.

- [ ] **Step 3: Run the full suite**

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```
Expected: PASS.

- [ ] **Step 4: Ask for approval, then commit (if any cleanup happened)**

```bash
git status
git add -u
git commit -m "clean up trailing ExtendedPhaseSpace references"
```

---

## Final verification checklist

- [ ] All tests pass: `julia --project=. -e 'using Pkg; Pkg.test()'`.
- [ ] `git grep -n 'ExtendedPhaseSpace\|minimize_simple_geometric_action\|proper_sgMAM_system' -- ":!docs/superpowers/"` returns empty.
- [ ] `FreidlinWentzellHamiltonian` is exported; `ExtendedPhaseSpace` is not defined.
- [ ] `minimize_geometric_action` accepts `::CoupledSDEs` and `::FreidlinWentzellHamiltonian` only; throws `MethodError` on `::CoupledODEs`.
- [ ] Rank-deficient SDE construction throws an `ArgumentError` containing the string `"rank-deficient"`.
- [ ] State-dependent diagonal and general multiplicative converge in both gMAM and sgMAM.
- [ ] `fw_action`, `geometric_action` give numerical agreement with analytic baselines on the 1D OU multiplicative test.
- [ ] `@inferred` passes on `update_p!`, `update_x!`, `geometric_gradient_step!`.
- [ ] Allocation regression tests pass.
- [ ] Docs build cleanly.
- [ ] Benchmarks load and run.
