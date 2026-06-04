# Rates, distributions & the generator

This page covers the quantities you can compute from the **infinitesimal generator** of an SDE: stationary and quasi-stationary distributions, mean first-passage times, metastable timescales, and the time evolution of densities. All of them are linear-algebra problems on a sparse discretisation of the generator on a Cartesian grid (the engine, described in the second half of this page). It supports any state-space dimension and any [`CoupledSDEs`](@ref) with diagonal additive noise.

## What you can compute

Every observable below is solved on a built [`DiffusionGenerator`](@ref) (see [The generator](@ref "The generator (engine)")). Building the generator is the expensive step (one sweep of the grid); reuse it across analyses.

```julia
grid = CartesianGrid((-2.0, 2.0, 81), (-2.0, 2.0, 81))
sys  = CoupledSDEs((u, p, t) -> SA[u[2], -u[1] + u[1]^3 - 0.5*u[2]],
                   [0.0, 0.0]; noise_strength = 0.4)
gen  = DiffusionGenerator(sys, grid)
```

The `A`, `B`, `target` arguments below accept a predicate `x -> Bool`, a `BitVector` of length `ncells(grid)`, or a vector of linear cell indices; the predicate builders [`ball`, `cuboid`, `sublevel`](@ref "Predicate and reshape helpers") construct common shapes.

By default the analyses use a sparse direct LU solve. For large grids where direct factorisation runs out of memory, pass an `alg` keyword with any [LinearSolve.jl](https://docs.sciml.ai/LinearSolve/stable/) algorithm, typically a Krylov method:

```julia
using LinearSolve: KrylovJL_GMRES
ρ = stationary_distribution(gen, KrylovJL_GMRES())
τ = mean_first_passage_time(gen, x -> x[1] > 1; alg = KrylovJL_GMRES())
```

### Stationary and quasi-stationary distributions

The **invariant density** ``\rho`` is the stationary solution of the Fokker-Planck equation (``\mathcal{L}^{*}\rho = 0`` with ``\int\rho = 1``); it is the long-time distribution of the system. The **quasi-stationary distribution** (QSD) is the analogous long-lived distribution *conditioned on not having left* a chosen metastable basin, obtained from the leading eigenmode of the generator restricted to that basin with absorbing boundary.

```julia
ρ   = stationary_distribution(gen)
qsd = quasi_stationary_distribution(gen, x -> x[1] < 0)   # QSD of the left basin
```

```@docs; canonical=false
stationary_distribution
quasi_stationary_distribution
```

### Mean first-passage time

The mean first-passage time ``\tau(x)`` to a target set is the expected time until a trajectory started at ``x`` first hits the target. Its variance is also available.

```julia
using CriticalTransitions: ball, reshape_to_grid

τ = mean_first_passage_time(gen, x -> x[1] > 1)
τ = mean_first_passage_time(gen, ball((1.0, 0.0), 0.25))

τ_grid = reshape_to_grid(τ, gen)            # vector → grid shape for plotting
heatmap(grid.centers[1], grid.centers[2], τ_grid')
```

```@docs; canonical=false
mean_first_passage_time
first_passage_variance
```

### Metastable timescales

The slowest-decaying eigenmodes of the generator encode metastable timescales (one over the spectral gap is the mixing time); the corresponding eigenvectors are the canonical reaction coordinates / metastable-state indicators used in Markov state model decomposition.

```julia
λ, V = eigenmodes(gen, 5)        # 5 slowest eigenpairs
@show λ                          # λ[1] = 0; λ[2:end] ≤ 0
metastable_timescale = -1 / real(λ[2])
```

```@docs; canonical=false
eigenmodes
```

### Time evolution of densities

The discrete Fokker-Planck equation `dρ/dt = Qᵀ ρ` is propagated via the matrix-exponential action `ρ(t) = exp(t·Qᵀ)·ρ₀`. The signature mirrors `DynamicalSystemsBase.trajectory`: total time `T`, optional sampling interval `Δt`, optional transient `Ttr`.

```julia
using CriticalTransitions: ncells, cell_volume

ρ_0 = zeros(ncells(gen)); ρ_0[i_start] = 1 / cell_volume(gen)   # point mass

ρs, t = propagate_density(gen, 5.0, ρ_0)               # endpoints only: t = [0, 5]
ρs, t = propagate_density(gen, 5.0, ρ_0; Δt = 0.1)     # snapshots every 0.1
ρs, t = propagate_density(gen, 5.0, ρ_0; Δt = 0.1, Ttr = 2.0)   # skip transient
```

`ρs[:, i]` is the density at `t[i]`. Internal substepping is adaptive; the `tol`, `m`, and `adaptive` kwargs control the Krylov-subspace matrix-exponential algorithm.

```@docs; canonical=false
propagate_density
```

## The generator (engine)

### Theory primer

We work with an Itô diffusion in ``\mathbb{R}^D``,

```math
\mathrm{d}X_t = b(X_t)\,\mathrm{d}t + \sigma(X_t)\,\mathrm{d}W_t .
```

The **infinitesimal generator** of the diffusion is the second-order differential operator

```math
\mathcal{L}\,u = b\!\cdot\!\nabla u + \tfrac{1}{2}\,\mathrm{tr}(\sigma\sigma^{\!\top}\nabla^{\!2} u) .
```

It is the **backward Kolmogorov operator**, governing how expectations of observables evolve in time (``\partial_t \mathbb{E}[u(X_t)] = \mathcal{L}u``). Its ``L^2``-adjoint ``\mathcal{L}^{*}`` is the **Fokker-Planck operator**, governing the probability density (``\partial_t \rho = \mathcal{L}^{*}\rho``).

After discretisation, ``\mathcal{L}`` becomes a sparse matrix ``Q``; its transpose ``Q^\top`` discretises ``\mathcal{L}^{*}``. The same matrix wears different names depending on convention:

| Operator | Acts on | Accessor |
|---|---|---|
| backward Kolmogorov / generator / rate matrix (= ``Q``) | observables | [`rate_matrix`](@ref) |
| (negative) M-matrix (= ``-Q``) | observables, PDE form | [`m_matrix`](@ref) |
| forward Kolmogorov / Fokker-Planck (= ``Q^\top``) | probability densities | [`fokker_planck_operator`](@ref) |

`rate_matrix` and `m_matrix` differ only in sign convention; the Fokker-Planck operator is the matrix transpose, equivalent to the ``L^2``-adjoint of the generator (for the uniform cell-volume inner product). The observables in the first half of this page are linear boundary-value problems on ``\mathcal{L}`` or its adjoint: the invariant density solves ``Q^\top\rho = 0`` with normalisation, and the mean first-passage time to a target ``T`` solves ``\mathcal{L}\tau = -1`` on ``T^{\mathrm{c}}`` with ``\tau|_T = 0``.

### Discretisation

We approximate the diffusion by a continuous-time Markov chain on a uniform Cartesian grid: each cell becomes a state of the chain; the generator ``Q`` collects transition rates between adjacent cells.

The face rates come from the **Scharfetter-Gummel** exponential-fitting stencil. It interpolates between centered differences (small Péclet) and upwind (large Péclet), and **always preserves the M-matrix property**, which is what guarantees solutions of ``Qq = 0`` numerically stay in the right range even at large drift / small noise.

The outer faces of the grid can be treated three different ways, set via the `bc` keyword to [`DiffusionGenerator`](@ref):

| `bc` | Meaning |
|---|---|
| [`Reflecting`](@ref)`()` (default) | outer faces omitted, the chain bounces off (no-flux Neumann). |
| [`Periodic`](@ref)`()` | outer faces wrap to the opposite end of the same axis. Use for systems on a torus / ring (angle variables). |
| [`Absorbing`](@ref)`()` | boundary cells leak probability through the outer faces; `Q` becomes a sub-generator. Use for survival problems where trajectories that leave the grid are removed. |

A single instance applies to every axis; pass a `D`-tuple for per-axis control:

```julia
gen = DiffusionGenerator(sys, grid; bc = Periodic())
gen = DiffusionGenerator(sys, grid; bc = (Reflecting(), Periodic()))
```

The boundary condition is part of the generator's type (a `DiffusionGenerator{2, Tuple{Reflecting, Periodic}}` is a different concrete type from `DiffusionGenerator{2, Tuple{Periodic, Periodic}}`), which lets the compiler specialize the assembly and analysis code to the specific BC. The covariance must be **diagonal** (axis-aligned noise); off-diagonal entries raise an `ArgumentError` at construction.

```@docs; canonical=false
CriticalTransitions.BoundaryCondition
Reflecting
Periodic
Absorbing
```

### Building the generator

```julia
gen = DiffusionGenerator(sys, grid)
```

```@docs; canonical=false
DiffusionGenerator
rate_matrix
m_matrix
fokker_planck_operator
```

### The grid

Each axis is `(lo, hi, N)`; per-axis spacing may differ.

```julia
grid = CartesianGrid((-2.0, 2.0, 81), (-2.0, 2.0, 81))   # 2D, 81×81 cells
```

```@docs; canonical=false
CartesianGrid
```

A handful of small accessors are available (non-exported, use them via `CriticalTransitions.ncells(grid)` or import explicitly):

```@docs; canonical=false
CriticalTransitions.ncells
CriticalTransitions.cell_volume
CriticalTransitions.cell_center
```

### Predicate and reshape helpers

Predicate builders for the `A` / `B` / `target` arguments and a vector to grid-shape reshaper, kept **non-exported** to keep the top-level namespace tight. Import explicitly:

```julia
using CriticalTransitions: ball, cuboid, sublevel, reshape_to_grid
```

```@docs; canonical=false
CriticalTransitions.ball
CriticalTransitions.cuboid
CriticalTransitions.sublevel
CriticalTransitions.reshape_to_grid
```

## Limitations

The implementation targets one specific point in design space: uniform Cartesian grid, Scharfetter-Gummel finite volume, sparse direct LU.

| | Status |
|---|---|
| State-space dimension | works for any ``D``, but practical only for ``D \leq 3`` (``\sim 10^6`` cells); marginal at ``D = 4``; infeasible at ``D \geq 5`` with uniform grids |
| Noise covariance | must be **diagonal** |
| Geometry | uniform Cartesian only; curved boundaries are staircased, no AMR |
| Boundary conditions | reflecting (default), periodic, or absorbing, per-axis |
| Time-dependence | autonomous only |
| Linear solver | sparse direct LU by default; opt-in Krylov via the `alg` kwarg (any LinearSolve.jl algorithm) |
