# Transition Path Theory

Transition Path Theory (TPT) studies **rare transitions** between metastable
states of an SDE [VandenEijndenTransition2006, VandenEijndenTransition2010, VandenEijndenTowards2006](@cite).
Given two disjoint sets ``A`` and ``B``, it tells you *how*, *where*, and
*how often* a noise-driven trajectory crosses from one to the other —
without simulating the (typically exponentially long) waiting times.

This package builds the discrete backward Kolmogorov operator on a uniform
Cartesian grid via the Scharfetter–Gummel finite-volume scheme, and solves
the resulting sparse linear systems. It works in any state-space dimension
and for any [`CoupledSDEs`](@ref) with diagonal additive noise.

---

## 1. Theory primer

### Setting

We work with an Itô diffusion in ``\mathbb{R}^D``,

```math
\mathrm{d}X_t = b(X_t)\,\mathrm{d}t + \sigma(X_t)\,\mathrm{d}W_t .
```

The drift has two metastable wells separated by a barrier. Pick disjoint
sets ``A`` and ``B`` enclosing them. TPT focuses on **reactive
trajectories**: pieces of the path that have just left ``A`` and will hit
``B`` *before* returning.

### The generator

The **infinitesimal generator** of the diffusion,

```math
\mathcal{L}\,u = b\!\cdot\!\nabla u + \tfrac{1}{2}\,\mathrm{tr}(\sigma\sigma^{\!\top}\nabla^{\!2} u),
```

discretises into a sparse matrix ``Q``. The same matrix wears different
names depending on convention:

| Name | Convention |
|---|---|
| **rate matrix / Q-matrix** | off-diagonals ``\geq 0``, row sums zero |
| **(negative) M-matrix** | the M-matrix is ``-Q`` |

Both views are exposed through [`rate_matrix`](@ref) and [`m_matrix`](@ref).

### Committors

The **forward committor** ``q^{+}(x)`` is the probability of reaching ``B``
before ``A`` from ``x``. It solves the boundary-value problem

```math
\mathcal{L}\,q^{+} = 0 \text{ on } (A\cup B)^{\mathrm{c}}, \qquad
q^{+}|_A = 0,\quad q^{+}|_B = 1 .
```

The **backward committor** ``q^{-}(x)`` is the probability that ``x`` last
visited ``A`` rather than ``B``, looking backwards in time. It solves the
same equation on the time-reversed process. For reversible systems
``q^{-} = 1 - q^{+}``; in general the two are independent.

### Reactive observables

Combining the committors with the **invariant density** ``\rho`` gives:

- **reactive density** ``\rho_R = \rho\,q^{+}\,q^{-}`` — probability
  density of being on a reactive piece;
- **reactive rate** ``k_{AB}`` — events per unit time, the flux of
  probability leaving ``A`` along reactive trajectories;
- **reactive current** ``J_R`` — vector field whose streamlines reveal
  the dominant transition channels.

### Mean first-passage time

A natural companion to the committor is ``\tau(x)``, the expected time to
first reach a target ``T``, solving ``\mathcal{L}\tau = -1`` on ``T^{\mathrm c}``
with ``\tau|_T = 0``. Useful well beyond TPT.

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

## 3. Discretisation

We approximate the diffusion by a continuous-time Markov chain on a
uniform Cartesian grid. Each cell becomes a state of the chain; the
generator ``Q`` collects transition rates between adjacent cells.

The face rates come from the **Scharfetter–Gummel** exponential-fitting
stencil. It interpolates between centered differences (small Péclet) and
upwind (large Péclet), and **always preserves the M-matrix property** —
which is what guarantees committors numerically stay in ``[0, 1]``, even
at large drift / small noise.

Outer-grid faces are simply omitted, giving **reflecting (no-flux Neumann)**
boundary conditions. Choose the grid so that ``A``, ``B``, and the
transition region all sit comfortably inside the box.

The covariance matrix must be **diagonal** (axis-aligned noise);
off-diagonal entries raise an `ArgumentError` at construction.

---

## 4. The API

The API is layered. Stop at any level and use the output as you wish, or
compose upward to the bundled [`ReactiveTransition`](@ref) result.

```text
Layer 4: ReactiveTransition  ← bundled TPT result
                  ▲
Layer 3: stationary_distribution, forward_committor, backward_committor,
         mean_first_passage_time
                  ▲
Layer 2: DiffusionGenerator  ← the discretised operator
                  ▲
Layer 1: CartesianGrid       ← the geometry
```

### Layer 1 — geometry

```julia
grid = CartesianGrid((-2.0, 2.0, 81), (-2.0, 2.0, 81))   # 2D, 81×81 cells
```

Each axis is `(lo, hi, N)`. Per-axis spacing may differ.

```@docs
CartesianGrid
```

### Layer 2 — the discretised operator

```julia
sys = CoupledSDEs((u, p, t) -> SA[u[2], -u[1] + u[1]^3 - 0.5*u[2]],
                  [0.0, 0.0]; noise_strength = 0.4)
gen = DiffusionGenerator(sys, grid)
```

Building `gen` is the expensive step (one sweep of the grid). Reuse it
across multiple analyses or sweeps over ``(A, B)``.

```@docs
DiffusionGenerator
rate_matrix
m_matrix
```

### Layer 3 — analyses

The `A`, `B`, `target` arguments accept a predicate `x -> Bool`, a
`BitVector` of length `ncells(grid)`, or a vector of linear cell indices.

```julia
ρ  = stationary_distribution(gen)
q⁺ = forward_committor(gen,  x -> x[1] < -1, x -> x[1] > 1)
q⁻ = backward_committor(gen, x -> x[1] < -1, x -> x[1] > 1)
τ  = mean_first_passage_time(gen, x -> x[1] > 1)
```

`backward_committor` defaults to the discrete adjoint of `gen.Q`. Pass
`reverse = gen_rev::DiffusionGenerator` to use a separately discretised
reverse-drift operator instead.

```@docs
stationary_distribution
forward_committor
backward_committor
mean_first_passage_time
```

### Layer 4 — bundled result

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

`forward_committor`, `backward_committor`, `stationary_distribution` also
dispatch on `ReactiveTransition` and return the cached fields.

```@docs
ReactiveTransition
reactive_rate
reactive_density
reactive_current
probability_reactive
probability_last_A
```

---

## 5. Complete example

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

## 6. Tips

- **Choose ``A`` and ``B`` generously** around the wells. Too small → noisy;
  overlapping the saddle → distorted current. A ball of radius
  ``\sim \sqrt{\sigma^2/\kappa}`` (with ``\kappa`` the well curvature) is
  a reasonable default.
- **The rate is in SDE time units.** Dividing by an appropriate occupation
  probability gives the mean transition time; see
  [VandenEijndenTransition2010](@cite) for the precise relation.
- **Sweeping ``(A, B)``?** Build [`DiffusionGenerator`](@ref) once and
  reuse it.

---

## 7. Limitations

The implementation targets one specific point in design space — uniform
Cartesian grid, Scharfetter–Gummel finite volume, sparse direct solver.

| | Status |
|---|---|
| State-space dimension | works for any ``D``, but practical only for ``D \leq 3`` (``\sim 10^6`` cells); marginal at ``D = 4``; infeasible at ``D \geq 5`` with uniform grids |
| Noise covariance | must be **diagonal** |
| Geometry | uniform Cartesian only — curved boundaries are staircased, no AMR |
| Boundary conditions | **reflecting** only (no absorbing or periodic) |
| Time-dependence | autonomous only |
| Linear solver | sparse direct LU; no Krylov |
| Small-noise regime | sharp committor layers near the saddle force globally fine grids — sgMAM / string method are the better tool there |

---

## References

```@bibliography
Pages = ["transition_path_theory.md"]
```
