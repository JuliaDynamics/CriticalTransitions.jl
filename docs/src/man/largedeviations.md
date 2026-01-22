# Large deviation theory

This section applies results of large deviation theory (LDT), particularly action minimization problems for computing most probable transition paths in stochastic dynamical systems driven by weak noise. For a description of the theory, see [freidlin_random_1998](@citet) and [borner_climate_2025](@cite). An overview of numerical methods applying LDT is given in [grafke_numerical_2019](@citet).

!!! info
    The methods in this section apply to ``D``-dimensional stochastic dynamical systems of the form
    ```math
    \text{d} \mathbf{x} = \mathbf{b} (\mathbf{x}) \text{d}t + \sigma \mathbf{\Sigma} \text{d}\mathbf{W}_t \,,
    ```
    where the drift field ``\mathbf{b}`` may be non-gradient but the noise term must consist of Gaussian noise (``\mathbf{W}_t`` is a ``D``-dimensional vector of independent standard Wiener processes) and a constant covariance matrix ``\mathbf{Q} = \mathbf{\Sigma}\mathbf{\Sigma}^\top``.

    This is a special case of the broader class of noise types supported by [`CoupledSDEs`](@ref).

## Action minimizers
Several methods have been proposed to calculate transition paths that minimize a given [action functional](@ref "Action functionals"). In the weak-noise limit, this minimum action path (or instanton) corresponds to the most probable transition path. While the minimum action method (MAM) is the most basic version, it is often beneficial to minimize the [geometric action](@ref "Geometric Freidlin-Wentzell action") via a time-independent version called gMAM. The problem can also be cast in a Hamiltonian form, implemented as simple gMAM (sgMAM), which can have numerical advantages.

These methods apply to non-gradient systems driven by Gaussian noise. In gradient systems, minimum action paths between attractors coincide with heteroclinic orbits, which can be computed via the so-called string method.

To summarize, the following methods are available:

- Minimum action method [(MAM)](@ref "Minimum action method (MAM)")
- Geometric minimum action method [(gMAM)](@ref "Geometric minimum action method (gMAM)")
- Simple geometric minimum action method [(sgMAM)](@ref "Simple geometric minimum action method (sgMAM)")
- [String method](@ref)

### Minimum action method (MAM)
Minimization of the specified action functional using the optimization algorithm of `Optimization.jl`. See also [e_minimum_2004](@citet).

```@docs
min_action_method
```

### Geometric minimum action method (gMAM)
Minimization of the geometric action following [heymann_pathways_2008, heymann_geometric_2008](@citet). gMAM reformulates MAM to avoid double optimization of both the action and the transition time. It achieves this by using a [geometric action](@ref "Geometric Freidlin-Wentzell action") functional that is independent of the time parametrization of the path. This reparameterization invariance makes the method more robust and computationally efficient, particularly for multiscale systems.

```@docs
geometric_min_action_method
```

### Simple geometric minimum action method (sgMAM)
Simplified minimization of the geometric action following [grafke_long_2017](@citet).
The simple gMAM reduces the complexity of the original gMAM by requiring only first-order derivatives of the underlying Hamiltonian optimization formulation. This simplifies the numerical treatment and computational complexity.

The implementation below performs a constrained gradient descent assuming an autonomous system with additive Gaussian noise.
```@docs
simple_geometric_min_action_method
ExtendedPhaseSpace
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
## Performance notes (sgMAM)

This page documents a common performance pitfall when using the simplified geometric minimum action method (sgMAM) with a general-purpose [`CoupledSDEs`](@ref) versus using a *hardcoded* extended phase space system.

### Summary

- sgMAM evolves a discretized path using Hamiltonian derivatives `H_p(x, p)` and `H_x(x, p)` in an extended phase space.
- If you construct the extended system from a `CoupledSDEs` via `ExtendedPhaseSpace(ds)`, the default implementation typically evaluates the drift and Jacobian *pointwise* along the path and can allocate heavily.
- A hardcoded `ExtendedPhaseSpace(H_x, H_p)` that uses analytic expressions (and operates on the full `D×Nt` path matrix) is usually much faster and much more allocation-friendly.

This often shows up as a large gap in allocations (and wall time) for the same `Nt` and number of iterations.

### Why `ExtendedPhaseSpace(ds::CoupledSDEs)` can allocate more

#### 1) Generic construction relies on the system Jacobian

When you call `simple_geometric_min_action_method(ds::ContinuousTimeDynamicalSystem, ...)`, the method internally constructs an [`ExtendedPhaseSpace`](@ref) from `ds`.
For `ds::CoupledSDEs`, this route uses `dynamic_rule(ds)` and `jacobian(ds)`.

If you did not provide an analytic Jacobian, `jacobian(ds)` is typically produced via automatic differentiation, and its evaluation may allocate (and is inherently more expensive than a hardcoded formula).

#### 2) The default `H_x` implementation can allocate inside tight loops

The convenience implementation of `ExtendedPhaseSpace(ds)` computes

\[ H_x(x, p) = J(x)^\top p \]

by forming a Jacobian matrix `J(x)` for each path point and then taking dot products.
If the Jacobian columns are accessed as `jax[:, idc]` without views, that column extraction creates temporary vectors repeatedly.

This is great for correctness and ease of use, but it is not tuned for allocation-free inner loops.

### Why a hardcoded extended phase space system is faster

For sgMAM you can often write `H_p`/`H_x` to operate on a full `D×Nt` matrix at once, using analytic expressions and broadcast fusion.

This tends to:

- Reduce Julia-level dispatch overhead (fewer calls into user code).
- Avoid per-point temporary vectors.
- Let broadcast fusion and BLAS-friendly operations kick in.

A particularly fast pattern is a “hardcoded” `H_p(x, p)` / `H_x(x, p)` pair that:

- accepts `x` and `p` as matrices (`D×Nt`),
- uses `eachrow(x)` / `eachrow(p)` and broadcasted formulas,
- returns a `D×Nt` matrix.

### A representative benchmark

The following benchmark structure (adapted from JuliaDynamics/CriticalTransitions.jl issue #145) compares:

- a hardcoded `ExtendedPhaseSpace` system used with [`simple_geometric_min_action_method`](@ref)
- an `ExtendedPhaseSpace(ds)` constructed from a `CoupledSDEs` drift (and typically an AD Jacobian)

```julia
using CriticalTransitions
using BenchmarkTools

# ... define fu, fv, df* and H_x/H_p as in the issue ...

sys = ExtendedPhaseSpace{false,2}(H_x, H_p)

function KPO(x, p, t)
    u, v = x
    return [fu(u, v), fv(u, v)]
end

ds = CoupledSDEs(KPO, zeros(2), ())

# "generic" extended phase space from CoupledSDEs
sys_generic = ExtendedPhaseSpace(ds)

# initial path matrix x_initial: size (D, Nt)

@btime simple_geometric_min_action_method($sys,         $x_initial; iterations=100, ϵ=0.5, show_progress=false)
@btime simple_geometric_min_action_method($sys_generic, $x_initial; iterations=100, ϵ=0.5, show_progress=false)
```

On one machine the hardcoded version was noticeably faster and allocated much less.

Exact numbers depend on Julia version, `Nt`, CPU, and how the drift is written, but the allocation gap is the key signal.

### Practical recommendations

- For sgMAM performance, prefer a hardcoded [`ExtendedPhaseSpace`](@ref) with analytic, matrix-based `H_p`/`H_x`.
- If you start from `CoupledSDEs`, consider providing an analytic Jacobian (when possible) so `ExtendedPhaseSpace(ds)` does not fall back to an expensive AD Jacobian.
- Use `@btime` to confirm allocations; big allocation counts are usually the main reason for slowdowns here.

As a brief aside: the same general principle (vectorized, allocation-free evaluation) also tends to make [`string_method`](@ref) faster when used with `ExtendedPhaseSpace`.
