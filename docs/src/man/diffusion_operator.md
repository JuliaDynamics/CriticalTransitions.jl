# Diffusion operator

This page documents how CriticalTransitions discretises an SDE on a
Cartesian grid into a sparse linear operator, and the analyses you can
compute on that operator: invariant density, mean first-passage time,
spectral decomposition, density propagation. It supports any state-space
dimension and any [`CoupledSDEs`](@ref) with diagonal additive noise.

---

## Theory primer

### Setting

We work with an Itô diffusion in ``\mathbb{R}^D``,

```math
\mathrm{d}X_t = b(X_t)\,\mathrm{d}t + \sigma(X_t)\,\mathrm{d}W_t .
```

### The generator

The **infinitesimal generator** of the diffusion is the second-order
differential operator

```math
\mathcal{L}\,u = b\!\cdot\!\nabla u + \tfrac{1}{2}\,\mathrm{tr}(\sigma\sigma^{\!\top}\nabla^{\!2} u) .
```

It is the **backward Kolmogorov operator**, governing how expectations of
observables evolve in time (``\partial_t \mathbb{E}[u(X_t)] = \mathcal{L}u``).
Its ``L^2``-adjoint ``\mathcal{L}^{*}`` is the **Fokker–Planck operator**,
governing the probability density (``\partial_t \rho = \mathcal{L}^{*}\rho``).

After discretisation, ``\mathcal{L}`` becomes a sparse matrix ``Q``; its
transpose ``Q^\top`` discretises ``\mathcal{L}^{*}``. The same matrix
wears different names depending on convention:

| Operator | Acts on | Accessor |
|---|---|---|
| backward Kolmogorov / generator / rate matrix (= ``Q``) | observables | [`rate_matrix`](@ref) |
| (negative) M-matrix (= ``-Q``) | observables, PDE form | [`m_matrix`](@ref) |
| forward Kolmogorov / Fokker–Planck (= ``Q^\top``) | probability densities | [`fokker_planck_operator`](@ref) |

`rate_matrix` and `m_matrix` differ only in sign convention; the
Fokker–Planck operator is the matrix transpose, equivalent to the
``L^2``-adjoint of the generator (for the uniform cell-volume inner
product).

### Boundary-value problems on the generator

A wide class of useful quantities solve linear BVPs on ``\mathcal{L}`` or
its adjoint. The ones implemented here:

- **Invariant density** ``\rho``: stationary solution of
  ``\mathcal{L}^{*}\rho = 0`` with ``\int \rho = 1``. Discretely:
  ``Q^\top \rho = 0`` with normalisation.
- **Mean first-passage time** ``\tau(x)`` to a target ``T``: expected time
  until ``X_t`` first hits ``T``. Solves ``\mathcal{L}\tau = -1`` on
  ``T^{\mathrm{c}}`` with ``\tau|_T = 0``.

---

## Discretisation

We approximate the diffusion by a continuous-time Markov chain on a
uniform Cartesian grid: each cell becomes a state of the chain; the
generator ``Q`` collects transition rates between adjacent cells.

The face rates come from the **Scharfetter–Gummel** exponential-fitting
stencil. It interpolates between centered differences (small Péclet) and
upwind (large Péclet), and **always preserves the M-matrix property** —
which is what guarantees solutions of ``Qq = 0`` numerically stay in the
right range, even at large drift / small noise.

The outer faces of the grid can be treated three different ways, set via
the `bc` keyword to [`DiffusionGenerator`](@ref):

| `bc` | Meaning |
|---|---|
| [`Reflecting`](@ref)`()` (default) | outer faces omitted — chain bounces off (no-flux Neumann). |
| [`Periodic`](@ref)`()` | outer faces wrap to the opposite end of the same axis. Use for systems on a torus / ring (angle variables). |
| [`Absorbing`](@ref)`()` | boundary cells leak probability through the outer faces; `Q` becomes a sub-generator. Use for survival problems where trajectories that leave the grid are removed. |

A single instance applies to every axis; pass a `D`-tuple for per-axis
control:

```julia
gen = DiffusionGenerator(sys, grid; bc = Periodic())
gen = DiffusionGenerator(sys, grid; bc = (Reflecting(), Periodic()))
```

The boundary condition is part of the generator's type — a
`DiffusionGenerator{2, Tuple{Reflecting, Periodic}}` is a different
concrete type from `DiffusionGenerator{2, Tuple{Periodic, Periodic}}` —
which lets the compiler specialize the assembly and analysis code to the
specific BC.

```@docs
CriticalTransitions.BoundaryCondition
Reflecting
Periodic
Absorbing
```

The covariance must be **diagonal** (axis-aligned noise); off-diagonal
entries raise an `ArgumentError` at construction.

---

## The API

The API is layered. Stop at any level and use the output as you wish.

```text
Layer 3: stationary_distribution, mean_first_passage_time,
         first_passage_variance, eigenmodes, propagate_density
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

A handful of small accessors are available (non-exported — use them via
`CriticalTransitions.ncells(grid)` or import explicitly):

```@docs
CriticalTransitions.ncells
CriticalTransitions.cell_volume
CriticalTransitions.cell_center
```

### Layer 2 — the discretised operator

```julia
sys = CoupledSDEs((u, p, t) -> SA[u[2], -u[1] + u[1]^3 - 0.5*u[2]],
                  [0.0, 0.0]; noise_strength = 0.4)
gen = DiffusionGenerator(sys, grid)
```

Building `gen` is the expensive step (one sweep of the grid). Reuse it
across multiple analyses or sweeps.

```@docs
DiffusionGenerator
rate_matrix
m_matrix
fokker_planck_operator
```

### Layer 3 — analyses

The `A`, `B`, `target` arguments accept a predicate `x -> Bool`, a
`BitVector` of length `ncells(grid)`, or a vector of linear cell indices.

```julia
ρ = stationary_distribution(gen)
τ = mean_first_passage_time(gen, x -> x[1] > 1)
```

The analyses default to a sparse direct LU solve. For large grids where
direct factorisation runs out of memory, pass an `alg` keyword with any
[LinearSolve.jl](https://docs.sciml.ai/LinearSolve/stable/) algorithm —
typically a Krylov method:

```julia
using LinearSolve: KrylovJL_GMRES
ρ = stationary_distribution(gen, KrylovJL_GMRES())
τ = mean_first_passage_time(gen, x -> x[1] > 1; alg = KrylovJL_GMRES())
```

```@docs
stationary_distribution
mean_first_passage_time
first_passage_variance
```

### Time evolution

The discrete Fokker–Planck equation `dρ/dt = Qᵀ ρ` is propagated via the
matrix-exponential action `ρ(t) = exp(t·Qᵀ)·ρ₀`. The signature mirrors
`DynamicalSystemsBase.trajectory` — total time `T`, optional sampling
interval `Δt`, optional transient `Ttr`:

```julia
using CriticalTransitions: ncells, cell_volume

ρ_0 = zeros(ncells(gen)); ρ_0[i_start] = 1 / cell_volume(gen)   # point mass

ρs, t = propagate_density(gen, 5.0, ρ_0)               # endpoints only: t = [0, 5]
ρs, t = propagate_density(gen, 5.0, ρ_0; Δt = 0.1)     # snapshots every 0.1
ρs, t = propagate_density(gen, 5.0, ρ_0; Δt = 0.1, Ttr = 2.0)   # skip transient
```

`ρs[:, i]` is the density at `t[i]`. Internal substepping is adaptive;
`tol`, `m`, and `adaptive` kwargs control the Krylov-subspace
matrix-exponential algorithm.

```@docs
propagate_density
```

### Spectral analysis

The slowest-decaying eigenmodes of `Q` encode metastable timescales (one
over the spectral gap is the mixing time) and the corresponding
eigenvectors are the canonical reaction coordinates / metastable-state
indicators used in Markov state model decomposition.

```julia
λ, V = eigenmodes(gen, 5)        # 5 slowest eigenpairs
@show λ                          # λ[1] = 0; λ[2:end] ≤ 0
metastable_timescale = -1 / real(λ[2])
```

```@docs
eigenmodes
```

---

## Convenience helpers

Predicate builders for `A` / `B` / `target` arguments and a vector →
grid-shape reshaper, kept **non-exported** to keep the top-level
namespace tight. Import explicitly:

```julia
using CriticalTransitions: ball, cuboid, sublevel, reshape_to_grid

target = ball((1.0, 0.0), 0.25)
τ = mean_first_passage_time(gen, target)

τ_grid = reshape_to_grid(τ, gen)
heatmap(grid.centers[1], grid.centers[2], τ_grid')
```

```@docs
CriticalTransitions.ball
CriticalTransitions.cuboid
CriticalTransitions.sublevel
CriticalTransitions.reshape_to_grid
```

---

## Limitations

The implementation targets one specific point in design space — uniform
Cartesian grid, Scharfetter–Gummel finite volume, sparse direct LU.

| | Status |
|---|---|
| State-space dimension | works for any ``D``, but practical only for ``D \leq 3`` (``\sim 10^6`` cells); marginal at ``D = 4``; infeasible at ``D \geq 5`` with uniform grids |
| Noise covariance | must be **diagonal** |
| Geometry | uniform Cartesian only — curved boundaries are staircased, no AMR |
| Boundary conditions | reflecting (default), periodic, or absorbing — per-axis |
| Time-dependence | autonomous only |
| Linear solver | sparse direct LU by default; opt-in Krylov via the `alg` kwarg (any LinearSolve.jl algorithm) |
