# OLIM Quasipotential Solver - Design

**Date:** 2026-05-13
**Target package:** `CriticalTransitions.jl`
**Status:** design, pre-implementation

## 1. Goal

Add a grid-based solver for the Freidlin-Wentzell quasipotential field $U(x)$ to `CriticalTransitions.jl`, using the **Ordered Line Integral Method** (OLIM) of Dahiya & Cameron (2018, *J. Sci. Comput.*, arXiv:1808.00562). The existing package only provides path-based methods (gMAM, sgMAM, string method) which compute one instanton at a time; OLIM produces the entire quasipotential field on a Cartesian domain in a single sweep.

The implementation must be dimension-generic via Julia's type system. It will be most accurate and most useful in $D = 2, 3$. For $D \geq 4$ a warning is emitted; for $D \geq 5$ the grid resolution per axis collapses fast enough that the result is informational only.

## 2. Scope

### In scope (v1)

- Single-attractor OLIM that fills the entire grid with $U_A(x)$ (the action from the supplied attractor).
- Compile-time dispatch on $D$ via `CoupledSDEs{IIP, D, I, P}` type parameter.
- Two quadrature rules: midpoint (`:midpoint`, 2nd order) and Simpson (`:simpson`, 4th order). Default `:simpson`.
- Simplex search on the accepted front for $D = 2$ (edge minimization) and $D = 3$ (triangle minimization). Vertex-only fallback for $D \geq 4$ via method dispatch.
- Analytic quadratic seeding via the continuous-time algebraic Riccati equation (CARE) in a small box around the attractor.
- `back_pointer` field on the result for instanton reconstruction by graph walk.
- Worked example (Maier-Stein 2D) under `examples/`.

### Out of scope (v1, planned follow-ups)

- Multi-attractor basin labeling and per-basin restriction.
- Saddle-stitched global quasipotential.
- Explicit caustic / Maxwell-set detection.
- Non-Cartesian or adaptive grids.
- Tensor-train / sparse-grid backends for higher dimensions.
- Right-endpoint quadrature (OUM-equivalent baseline) and other Cameron variants beyond MID/SIM.

## 3. Public API

```julia
struct QuasiPotential{D}
    U::Array{Float64, D}                       # +Inf only on unreached cells
    back_pointer::Array{CartesianIndex{D}, D}  # zero index for source / unreached
    source::CartesianIndex{D}
    grid::CartesianGrid{D}
end

"""
    olim(sys, grid, attractor; kwargs...) -> QuasiPotential{D}

Compute the Freidlin-Wentzell quasipotential field U_A(x) using the Ordered
Line Integral Method (Dahiya & Cameron 2018). The state dimension D is taken
from sys::CoupledSDEs{IIP, D, I, P} and must match grid::CartesianGrid{D}.

# Keyword arguments
- quadrature::Symbol = :simpson    -- :midpoint or :simpson
- band_radius::Int   = 2           -- K, the accepted-band radius in grid cells
- near_source_layers::Int = 3      -- size of the analytic CARE seed box; 0 disables
- verbose::Bool      = false
- show_progress::Bool = true

A warning is emitted for D > 4; the algorithm still runs but per-axis grid
resolution typically becomes too coarse to be useful past D = 4.
"""
function olim(sys::CoupledSDEs{IIP, D, I, P},
              grid::CartesianGrid{D},
              attractor::AbstractVector{<:Real};
              quadrature::Symbol = :simpson,
              band_radius::Int   = 2,
              near_source_layers::Int = 3,
              verbose::Bool      = false,
              show_progress::Bool = true) where {IIP, D, I, P}
    length(attractor) == D || throw(DimensionMismatch(
        "attractor has length $(length(attractor)) but sys has D=$D"))
    x_A = SVector{D, Float64}(attractor)
    return _olim_impl(sys, grid, x_A,
                      Val(quadrature), Val(band_radius), Val(near_source_layers),
                      verbose, show_progress)
end
```

Exports: `olim`, `QuasiPotential`.

The public function accepts any `AbstractVector{<:Real}` for ergonomics (users typically pass `[x, y]`), converts it once to `SVector{D, Float64}`, and wraps the `Symbol` / `Int` kwargs into `Val{quadrature}`, `Val{band_radius}`, `Val{near_source_layers}` before calling the internal type-stable `_olim_impl`. All inner functions specialize on those `Val` types and on the `SVector` element type.

## 4. Algorithm

OLIM is a Dijkstra-style one-pass sweep of the static Hamilton-Jacobi equation $H(x, \nabla U) = 0$ with $H(x, p) = p \cdot b(x) + \tfrac{1}{2} p \cdot a(x) p$, where $a = \sigma \sigma^\top$ is the noise covariance. Causality (the value at a Considered cell only depends on neighbors with smaller $U$) makes the one-pass sweep correct.

### 4.1 Outer loop

Standard Dijkstra sweep with three node statuses (`UNKNOWN`, `CONSIDERED`, `ACCEPTED`), a min-heap keyed by tentative $U$, and an accepted band of radius $K$ around each Considered point.

```
seed_near_source!(state, Val(D), Val(K_seed))            # analytic CARE box
push every CONSIDERED cell from the seeding step onto heap

while heap not empty:
    pop (_, lin) -> c
    status[c] = ACCEPTED

    for x in stencil(c, Val(K), bounds):                  # K-radius hyperbox around c
        status[x] == ACCEPTED && continue
        (newU, newback) = local_update(Val(D), x, state, L,
                                       Val(quad), Val(K))
        newU < U[x] || continue
        U[x]            = newU
        back_pointer[x] = newback
        if status[x] == UNKNOWN
            status[x] = CONSIDERED
            handles[lin(x)] = push!(heap, (newU, lin(x)))
        else            # already CONSIDERED
            update!(heap, handles[lin(x)], (newU, lin(x)))
        end
    end
end

return QuasiPotential(U, back_pointer, source, grid)
```

Two paths into the heap: `push!` for `UNKNOWN -> CONSIDERED` transitions, `update!` (decrease-key) for already-`CONSIDERED` cells whose tentative value improved. The handle map persists across both.

Heap: `DataStructures.MutableBinaryHeap{Tuple{Float64, Int}, FasterForward}` plus a `Dict{Int, Int}` from linear cell index to heap handle to support decrease-key.

Stencil: `CartesianIndices` over the fixed-size hyperbox of side $2K+1$ centered at the popped cell. With $K$ lifted to `Val{K}` the stencil is a compile-time-known offset list; for $K = 2$ in $D = 2$ this is 25 entries, fully unrollable.

### 4.2 Local update

The OLIM-specific part. For a Considered cell $x$:

1. Collect the set $\mathcal{A}(x)$ of `ACCEPTED` cells within Chebyshev distance $K$ of $x$.
2. Vertex candidates: for each $y \in \mathcal{A}(x)$, compute $\Phi_{\text{v}}(y) = U(y) + \int_0^1 L(y + s(x-y),\, x-y)\, ds$.
3. Simplex candidates (only when methods are defined, i.e. $D = 2$ or $3$): minimize $\Phi$ over straight segments from interior points of $(d{-}1)$-simplices on the accepted front.
4. Take the overall minimum.

Per-$D$ specialization is achieved by method dispatch on `Val{D}`:

```julia
add_simplex_candidates!(best, x, accepted, state, L, ::Val{quad}, ::Val{D}) where {quad, D} = best
add_simplex_candidates!(best, x, accepted, state, L, ::Val{quad}, ::Val{2}) where {quad} = ...
add_simplex_candidates!(best, x, accepted, state, L, ::Val{quad}, ::Val{3}) where {quad} = ...
```

The default method is a no-op (vertex-only path); specialized methods exist only for $D = 2$ and $D = 3$. The compiler picks the right one at each call site; there is no runtime branch on $D$.

### 4.3 Vertex candidate

Straight segment from accepted $y$ to Considered $x$, constant velocity $v = x - y$:

$$\Phi_{\text{v}}(y) = U(y) + \int_0^1 L\big(y + s v,\ v\big)\, ds,\qquad L(z, v) = \tfrac{1}{2}\big(v - b(z)\big)^\top a(z)^{-1} \big(v - b(z)\big).$$

### 4.4 Edge minimization (D = 2)

For an admissible pair $(y_0, y_1)$ of mutually adjacent accepted vertices, parametrize the back-point as $y(\lambda) = (1-\lambda) y_0 + \lambda y_1$ for $\lambda \in [0, 1]$. The objective is

$$\Phi(\lambda) = U_{\text{interp}}(\lambda) + \int_0^1 L\big(y(\lambda) + s(x - y(\lambda)),\ x - y(\lambda)\big)\, ds,$$

where $U_{\text{interp}}(\lambda) = (1-\lambda) U(y_0) + \lambda U(y_1)$ for `:midpoint`, or a quadratic in $\lambda$ using a third nearby accepted vertex for `:simpson`. Minimization on $[0, 1]$ via Brent's method, with $\lambda = 0$ and $\lambda = 1$ evaluated explicitly so corner minima are handled cleanly.

### 4.5 Triangle minimization (D = 3)

For an admissible triple $(y_0, y_1, y_2)$ of mutually adjacent accepted vertices, parametrize the back-point in barycentric coordinates:

$$y(\lambda_1, \lambda_2) = (1 - \lambda_1 - \lambda_2) y_0 + \lambda_1 y_1 + \lambda_2 y_2,\qquad \lambda_1, \lambda_2 \geq 0,\ \lambda_1 + \lambda_2 \leq 1.$$

2D constrained minimization on the unit triangle. Newton on the unconstrained gradient with projection back onto the simplex; fallback to a coarse grid scan plus Nelder-Mead if Newton stalls or leaves the feasible region. The three edges are always also evaluated as $D = 2$ problems (cheap, catches boundary minima).

### 4.6 Admissibility of simplices

A pair (or triple) is admissible iff:

- all vertices have `status == ACCEPTED`,
- they are mutually adjacent on the grid (Chebyshev distance 1 between every pair),
- they lie within Chebyshev radius $K$ of $x$,
- the segment / triangle does not contain $x$ (i.e. $x$ projects strictly outside on the "Unknown" side).

The admissibility predicate is also dispatched per $D$ via specialized iterators `front_edges`, `front_triangles`.

### 4.7 Back-pointer convention

When the winning candidate is a simplex (continuous interior minimizer), the back-pointer stored is the discrete vertex of that simplex closest in $\ell^2$ to the continuous minimizer. This keeps `back_pointer` a `CartesianIndex{D}` and lets users reconstruct a discrete MEP by chasing pointers. Users wanting sub-cell accuracy can refine the pointer-walked path with `minimize_geometric_action`.

## 5. Quadrature

Two rules, lifted to types via `Val{:midpoint}` and `Val{:simpson}`:

```julia
@inline line_integral(L::FWLagrangian{D}, y, v, ::Val{:midpoint}) where {D} =
    L(y + 0.5 * v, v)

@inline line_integral(L::FWLagrangian{D}, y, v, ::Val{:simpson}) where {D} =
    (L(y, v) + 4 * L(y + 0.5 * v, v) + L(y + v, v)) / 6
```

No symbol comparison at runtime; each method body inlines into the call site.

## 6. Lagrangian wrapper

```julia
struct FWLagrangian{D, B, A}
    b::B                # SVector{D, Float64} -> SVector{D, Float64}
    a_inv::A            # SMatrix{D,D,Float64,D*D} for constant noise; callable for state-dep
end

@inline function (L::FWLagrangian{D})(x::SVector{D, Float64},
                                       v::SVector{D, Float64}) where {D}
    Δ = v - L.b(x)
    return 0.5 * dot(Δ, _a_inv_at(L.a_inv, x), Δ)
end

@inline _a_inv_at(M::SMatrix, _) = M
@inline _a_inv_at(f::F, x) where {F<:Function} = f(x)
```

Construction extracts the drift via `dynamic_rule(sys)` and the diffusion via `covariance_matrix(sys)`. For an out-of-place `CoupledSDEs` the drift is used directly; for an in-place `CoupledSDEs` it is wrapped at construction in a small allocating closure that calls the in-place form into a stack-allocated `MVector{D, Float64}`. Either way the resulting `b::B` field of `FWLagrangian` has the OOP signature `SVector{D} -> SVector{D}`, so the `IIP` flag is consumed at construction time and does not need to live on the struct.

For constant noise the inverse is precomputed once as an `SMatrix{D,D,Float64,D*D}`; for state-dependent noise it is a closure over `sigma(x)`.

## 7. Near-source analytic seeding

Inside a box of side $2 K_{\text{seed}} + 1$ centered at the attractor cell, set $U(x) \approx (x - x_A)^\top Q (x - x_A)$ where $Q$ solves the continuous-time algebraic Riccati equation

$$A^\top Q + Q A + 2\, Q\, a\, Q = 0,\qquad A = \nabla b(x_A),\ a = \sigma \sigma^\top.$$

Solver: `MatrixEquations.arec` for $D \geq 3$; closed-form for $D = 2$. The drift Jacobian $A$ is obtained via `ForwardDiff.jacobian`.

Cells inside the box are marked `ACCEPTED` with their analytic $U$ values. The first sweep iteration picks up from the box edge.

Guard: if `eigvals(A)` has any non-negative real part, $x_A$ is not a stable fixed point of $b$; skip analytic seeding (emit `@info` when `verbose`), start the sweep from the source cell with $U = 0$ and `near_source_layers = 0` semantics.

`near_source_layers = 0` is supported and disables analytic seeding.

## 8. Compile-time dispatch summary

Every dimension- or quadrature-dependent code path is selected by Julia method dispatch:

| Decision | Dispatched on |
|---|---|
| state dimension | `Val{D}` (specialized methods for `Val{2}`, `Val{3}`, generic fallback) |
| quadrature rule | `Val{:midpoint}`, `Val{:simpson}` |
| band radius (stencil size) | `Val{K}` |
| near-source box size | `Val{K_{\text{seed}}}` |
| in-place vs out-of-place drift | `IIP` parameter inherited from `CoupledSDEs` |
| constant vs state-dependent noise | type of `FWLagrangian.a_inv` |

There is no `if D == ...`, `if quadrature == ...`, or `if band_radius == ...` branch at runtime inside any hot loop.

## 9. File layout

```
src/largedeviations/olim/
├── olim.jl              public API, QuasiPotential, _olim_impl
├── state.jl             OLIMState container (U, back_pointer, status, heap, handles)
├── sweep.jl             Dijkstra sweep loop
├── update.jl            local_update dispatcher, vertex_candidate
├── simplex.jl           edge_minimum (D=2), triangle_minimum (D=3), front_edges, front_triangles
├── quadrature.jl        line_integral methods
├── lagrangian.jl        FWLagrangian struct and constructors
└── source_seed.jl       Riccati solve, seed_near_source!
```

`src/CriticalTransitions.jl` adds `include("largedeviations/olim/olim.jl")` and exports `olim` and `QuasiPotential`. `olim/olim.jl` includes its siblings in dependency order.

## 10. Dependencies

New direct dependencies added to `Project.toml`:

| Package | Compat | Reason |
|---|---|---|
| `DataStructures` | ≥ 0.18 | `MutableBinaryHeap` for the Dijkstra priority queue with decrease-key |
| `MatrixEquations` | ≥ 2.0 | `arec` for the CARE-based near-source seed |
| `ForwardDiff` | already present | Jacobian of drift at the attractor |

`StaticArrays`, `LinearAlgebra`, `SparseArrays` are already direct deps.

## 11. Testing

`test/largedeviations/olim.jl`:

1. **Quadrature units** - integrate $L(s\,v, v)$ for analytic $b$, $a$; verify $O(h^2)$ for `:midpoint`, $O(h^4)$ for `:simpson`.
2. **Vertex update unit** - 2D quadratic well, analytic $U$; verify single-cell `vertex_candidate` to expected order.
3. **Edge minimization unit (D=2)** - synthetic problem with analytic optimal $\lambda$; Brent recovers it to `1e-10`.
4. **Triangle minimization unit (D=3)** - synthetic with known interior minimum.
5. **End-to-end gradient (2D double well)** - $U = 2V$ on a $100^2$ grid; relative error tolerance 5% with `:midpoint`, 1% with `:simpson`.
6. **End-to-end Maier-Stein (2D non-gradient)** - $U$ at the saddle vs the literature value (≈ 0.5); also vs `minimize_geometric_action`. Tolerance 2%.
7. **3D linear Ornstein-Uhlenbeck smoke** - analytic quadratic $U$; coarse $30^3$ grid; mid-grid relative error < 5%.
8. **D=5 warning** - tiny $15^5$ grid; `@test_logs (:warn, r"OLIM in D=5") olim(...)`.
9. **Back-pointer walk** - from the saddle cell, follow the `back_pointer` chain; verify termination at the source and monotone-decreasing $U$ along the trace.

Tests follow the existing `test/largedeviations/` style: `@testset` blocks, seeded RNG, `CoupledSDEs` constructors with out-of-place drifts returning `SVector`.

## 12. Example

`examples/olim_maierstein.jl`:

- Define the Maier-Stein system as a `CoupledSDEs` (matches `examples/gMAM_Maierstein.jl`).
- Build a $200^2$ `CartesianGrid` covering $[-1.5, 1.5]^2$.
- Run `olim(sys, grid, SVector(-1.0, 0.0))` from the left attractor.
- Visualize $U$ with `CairoMakie.contourf`, overlay the gMAM instanton, mark attractor and saddle.
- Print the OLIM action at the saddle cell and the gMAM action for comparison.

## 13. Documentation

- `docs/src/man/largedeviations.md` (or the analogous existing page) gains a section on OLIM with a one-screen usage block and a pointer to the Maier-Stein example.
- `docs/src/refs.bib` gains the Dahiya-Cameron 2018 entry and Cameron's earlier 2012 *Physica D* paper.
- `CHANGELOG.md` gets an entry under the next release.

## 14. Risks and known limitations

- **Memory** scales as $N^D$. The `D > 4` warning is informational; we do not refuse to run. Users on small machines are responsible for sizing the grid.
- **Heap operations** are $O(D \log N)$ per pop/push; the dominant cost in practice is the local update inside the stencil, not the heap.
- **Triangle minimization (D=3)** is the most fragile piece. The Newton-with-projection fallback is documented and tested on synthetic problems; expect to revisit if benchmarks reveal pathological cases.
- **State-dependent noise** is supported by the type system but not benchmarked beyond a smoke test in v1.
- **Caustics** are not detected. The sweep returns the viscosity solution (correct minimum), but users wanting to visualize Maxwell sets will need to layer that on themselves in a follow-up.

## 15. References

1. S. Dahiya and M. K. Cameron. *Ordered Line Integral Methods for Computing the Quasi-Potential.* Journal of Scientific Computing 75 (2018). arXiv:1808.00562. Cameron's accompanying code: https://math.umd.edu/~mariakc/OLIM.html
2. M. K. Cameron. *Finding the quasipotential for nongradient SDEs.* Physica D 241 (2012).
3. J. A. Sethian and A. Vladimirsky. *Ordered upwind methods for static Hamilton-Jacobi equations.* SIAM Journal on Numerical Analysis 41 (2003).
4. T. Grafke, T. Schäfer, and E. Vanden-Eijnden. *Long term effects of small random perturbations on dynamical systems: theoretical and computational tools.* In *Recent Progress and Modern Challenges in Applied Mathematics* (2019); arXiv:1604.03818.
