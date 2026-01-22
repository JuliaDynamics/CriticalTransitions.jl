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

| Method | Use when | Requirements (this package) | Not suitable when |
|---|---|---|---|
| **MAM** | You want a **minimum action path** for a specified travel time $T$ (FW/OM action minimization). | `CoupledSDEs` with **additive**, **invertible**, **autonomous** noise (constant covariance); discretized path uses **equispaced time**. | **Multiplicative/state-dependent noise**, **degenerate/non-invertible** noise, **non-autonomous** noise; when the transition time is unknown and you want a time-reparametrization-invariant formulation (prefer gMAM/sgMAM). |
| **gMAM** | You want a **time-reparametrization-invariant** minimum action path (no explicit optimization over $T$). | Same as MAM for `CoupledSDEs`: **additive**, **invertible**, **autonomous** noise (constant covariance). | **Multiplicative/state-dependent**, **degenerate**, or **non-autonomous** noise; if you need the Onsager--Machlup functional (only FW has a geometric formulation). |
| **sgMAM** | You want a **Hamiltonian/simple gMAM** formulation that can be efficient in practice. | Same as MAM **and** **diagonal** noise covariance (implementation restriction); assumes **additive** noise. | **Non-diagonal** covariance; **multiplicative/state-dependent** or **degenerate** noise; models where the required derivatives/Jacobian are not available or are too expensive. |
| **String method** | You want a **minimum energy path / heteroclinic orbit** driven by the deterministic drift (typical use: **gradient** systems). | Deterministic drift field (works for `ContinuousTimeDynamicalSystem`; does not rely on an SDE noise model). | In **non-gradient** systems if you need the *most probable* noise-induced transition path (string gives the deterministic heteroclinic orbit, which generally differs from the instanton). |

#### Variants and extensions
The literature contains a number of extensions of MAM-type methods that may be relevant depending on the model class and numerical difficulties. These variants are not currently implemented in `CriticalTransitions.jl`, but serve as useful pointers:

- **tMAM / optimal linear time scaling**: avoids explicit optimization over the transition time by introducing an optimal linear time scaling; can be combined with adaptivity in time discretization [wan_tmam_2015](@citet).
- **Adaptive MAM**: uses a moving-mesh strategy to concentrate grid points in dynamically important portions of the path, improving efficiency and robustness [zhou_adaptive_mam_2008](@citet).
- **Non-Gaussian (jump / Lévy) noise**: for systems driven by jump noise, the rate function and path optimization problem differ from the Freidlin--Wentzell diffusive setting; see e.g. an optimal-control-based approach in [wei_most_likely_jumps_2023](@citet).
- **Multiplicative / state-dependent noise**: extensions of geometric action minimization to degenerate or multiplicative noise are discussed in [grafke_small_random_2017](@citet); the current `sgMAM` implementation assumes additive noise.

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
