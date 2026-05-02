# Transition Path Theory

Transition Path Theory (TPT) studies **rare transitions** between
metastable states of an SDE
[VandenEijndenTransition2006, VandenEijndenTransition2010, VandenEijndenTowards2006](@cite).
Given two disjoint sets ``A`` and ``B``, it tells you *how*, *where*, and
*how often* a noise-driven trajectory crosses from one to the other —
without simulating the (typically exponentially long) waiting times.

TPT is built on top of the discretised generator and the standard
boundary-value problems on it (committors, invariant density). Read the
[diffusion operator manual](@ref Diffusion-operator) first if you have not
already.

---

## 1. Reactive trajectories and observables

The key objects are the **reactive trajectories**: pieces of a long path
that have just left ``A`` and will hit ``B`` *before* returning. TPT
characterises this sub-ensemble through three observables, all built from
the forward committor ``q^{+}``, the backward committor ``q^{-}``, and the
invariant density ``\rho``:

- **Reactive density** ``\rho_R = \rho\,q^{+}\,q^{-}`` — probability
  density of being on a reactive piece at a given time.
- **Reactive transition rate** ``k_{AB}`` — events per unit time, the flux
  of probability leaving ``A`` along reactive trajectories.
- **Reactive current** ``J_R`` — vector field whose streamlines reveal
  the dominant transition channels.

For reversible systems ``q^{-} = 1 - q^{+}``; in general the two are
independent and both must be computed.

The reactive current splits into a **reversible** (gradient-flow) part
and an **irreversible** (cyclic) part — the latter vanishes for systems
satisfying detailed balance and exposes the circulating component of the
flow for non-equilibrium systems. The pieces are accessed via
[`reactive_current_reversible`](@ref) and
[`reactive_current_irreversible`](@ref); they sum to the full
[`reactive_current`](@ref).

---

## 2. When TPT is the right tool

**Use TPT when** you have a time-homogeneous SDE with metastable wells, and
you want **spatial information** (where reactive paths concentrate, the
committor surface, transition channels) — not just a scalar rate.
Practical sweet spot: dimensionless barrier ``\Delta U / \sigma^2`` between
``\sim 1`` and ``10``, state-space dimension ``\lesssim 3``.

**Look elsewhere when:**

- **Very high barriers** (``\Delta U/\sigma^2 \gg 10``): use sgMAM / string
  method ([`minimize_action`](@ref), [`string_method`](@ref)) — they
  compute the most likely path and its rate directly.
- **Very shallow barriers**: just simulate ([`transitions`](@ref)).
- **Only need a rate**: use [`mean_first_passage_time`](@ref) or Kramers
  asymptotics.
- **Not a diffusion** (jumps, fractional noise, deterministic chaos):
  outside the SDE framework on which TPT is built.
- **High-dimensional** (``D \gtrsim 4``) with no reaction coordinate: a
  uniform grid is infeasible. The standard remedy is to first learn a 1–3D
  collective variable and *then* run TPT in CV space — outside this package.
- **Non-autonomous drift**: not implemented here.

A useful sanity check is to compare the TPT rate against an independent
estimate (Kramers, simulation, or large-deviation action).

---

## 3. The API

[`ReactiveTransition`](@ref) runs the full pipeline (``\rho``, both
committors, the rate) in one call and caches the results.

```julia
result = ReactiveTransition(sys, grid, A, B)         # or
result = ReactiveTransition(gen, A, B)               # reuse a generator
result = ReactiveTransition(sys, grid, A, B; reverse = sys_rev)

reactive_rate(result)                # cached
reactive_density(result)
J_nodes, J_faces = reactive_current(result)
probability_reactive(result)
probability_last_A(result)
```

`forward_committor`, `backward_committor`, and `stationary_distribution`
also dispatch on a `ReactiveTransition` and return the cached fields.
Predicates / `BitVector` / linear-index masks are accepted for `A` and `B`,
just as on a [`DiffusionGenerator`](@ref).

```@docs
ReactiveTransition
reactive_rate
reactive_density
reactive_current
reactive_current_reversible
reactive_current_irreversible
probability_reactive
probability_last_A
```

---

## 4. Complete example

A 2D Maier–Stein-like system with stable wells at ``(\pm 1, 0)``:

```julia
using CriticalTransitions
using LinearAlgebra: norm

b(u, p, t) = SA[u[1] - u[1]^3 - 10*u[1]*u[2]^2, -(1 + u[1]^2)*u[2]]
sys = CoupledSDEs(b, [0.0, 0.0]; noise_strength = 0.3)

grid = CartesianGrid((-1.6, 1.6, 121), (-1.0, 1.0, 81))
A = x -> norm(x .- SA[-1.0, 0.0]) < 0.25
B = x -> norm(x .- SA[ 1.0, 0.0]) < 0.25

result = ReactiveTransition(sys, grid, A, B)
@show reactive_rate(result)
J_nodes, _ = reactive_current(result)
```

To **sweep** over ``(A, B)`` choices without re-discretising:

```julia
gen = DiffusionGenerator(sys, grid)
for r in (0.2, 0.3, 0.4)
    A = x -> norm(x .- SA[-1.0, 0.0]) < r
    B = x -> norm(x .- SA[ 1.0, 0.0]) < r
    println("r = $r → k_AB = ", reactive_rate(ReactiveTransition(gen, A, B)))
end
```

A fully worked, plotted version is in the
[double-well example](@ref TPT_example).

---

## 5. Tips

- **Choose ``A`` and ``B`` generously** around the wells. Too small → noisy;
  overlapping the saddle → distorted current. A ball of radius
  ``\sim \sqrt{\sigma^2/\kappa}`` (with ``\kappa`` the well curvature) is
  a reasonable default.
- **The rate is in SDE time units.** Dividing by an appropriate occupation
  probability gives the mean transition time; see
  [VandenEijndenTransition2010](@cite) for the precise relation.
- **Sweeping ``(A, B)``?** Build [`DiffusionGenerator`](@ref) once and
  reuse it.
- **Small-noise regime.** As ``\sigma \to 0`` the committor develops a
  thin transition layer near the saddle — uniform grids must be globally
  fine, which is wasteful. sgMAM / string method are the better tool there.

---

## References

```@bibliography
Pages = ["transition_path_theory.md"]
```
