# OLIM Quasipotential Solver - Design

**Date:** 2026-05-13 (refreshed 2026-05-25 against `main` post #327, #329, #335, #338)
**Target package:** `CriticalTransitions.jl`
**Status:** design, pre-implementation

## 1. Goal

Add a grid-based solver for the Freidlin-Wentzell quasipotential field $U(x)$ to `CriticalTransitions.jl`, using the **Ordered Line Integral Method** (OLIM) of Dahiya & Cameron (2018, *J. Sci. Comput.*, arXiv:1808.00562). The existing package only provides path-based methods (gMAM, sgMAM, string method) which compute one instanton at a time; OLIM produces the entire quasipotential field on a Cartesian domain in a single sweep.

The implementation must be dimension-generic via Julia's type system. It will be most accurate and most useful in $D = 2, 3$. For $D \geq 4$ a warning is emitted; for $D \geq 5$ the grid resolution per axis collapses fast enough that the result is informational only.

## 2. Scope

### In scope (v1)

- Single-attractor OLIM that fills the entire grid with $U_A(x)$ (the action from the supplied attractor).
- Compile-time dispatch on $D$ via the `CoupledSDEs{IIP, D, I, P}` type parameter (DynamicalSystemsBase).
- Reuses the existing `CartesianGrid{D, T}` (from `src/diffusion_operator/cartesian_grid.jl`, introduced in #327) as the domain type.
- **Geometric** (arc-length / Maupertuis) Lagrangian $L_g(x, v) = \|v\|_Q \|b(x)\|_Q - \langle v, b(x)\rangle_Q$ with $Q = a/s$ trace-normalized via `_trace_normalized_a(sys)` to match `fw_action` conventions (`src/largedeviations/action.jl`, post-#329). Using $L_g$ (not the time-parametrized $\tfrac{1}{2}\|v-b\|_Q^2$) is what makes the one-pass Dijkstra sweep monotone and correct.
- Multiplicative (state-dependent) noise supported as a first-class case via `diffusion_function(sys)` (post-#338), with $Q(x)$ re-evaluated at every quadrature node.
- Simpson's rule (4th order) for the line integral along each candidate segment, paired with Hermite-cubic interpolation of $U$ along the front for consistent 4th-order accuracy.
- Four-status node tagging (`UNKNOWN`, `CONSIDERED`, `FRONT`, `ACCEPTED`) with an explicit `FRONT::BitArray{D}` so the 2-point edge / 3-point triangle update only iterates the frontier of accepted cells.
- Concrete simplex search on the front for $D = 2$ (edge minimization) and $D = 3$ (active-set KKT on the unit triangle). Vertex-only fallback for $D \geq 4$ via method dispatch.
- Analytic quadratic seeding via the continuous-time algebraic Riccati equation (CARE) in a small box around the attractor.
- Sub-cell `BackRef{D}` per cell (winning vertex pair + interior parameter `s`) for instanton reconstruction via graph walk.
- Worked example (Maier-Stein 2D) under `examples/`.

### Out of scope (v1, planned follow-ups)

- Multi-attractor basin labeling and per-basin restriction.
- Saddle-stitched global quasipotential.
- Explicit caustic / Maxwell-set detection.
- Non-Cartesian or adaptive grids.
- Tensor-train / sparse-grid backends for higher dimensions.
- Alternative quadrature rules (midpoint, right-endpoint / OUM-equivalent baseline) and other Cameron variants beyond Simpson.

## 3. Public API

```julia
"""
Sub-cell back-reference for instanton reconstruction.
For one-point updates (vertex), `v1 == v0` and `s == NaN32`.
For two-point (edge) updates, the predecessor lies at `(1-s)*v0 + s*v1`.
"""
struct BackRef{D}
    v0::CartesianIndex{D}
    v1::CartesianIndex{D}
    s::Float32
end

struct QuasiPotential{D, T <: AbstractFloat}
    U::Array{T, D}                       # +Inf only on unreached cells
    back_pointer::Array{BackRef{D}, D}   # zeroed BackRef for source / unreached
    source::CartesianIndex{D}
    grid::CartesianGrid{D, T}
end

"""
    quasipotential(sys, grid, attractor; kwargs...) -> QuasiPotential{D}

Compute the Freidlin-Wentzell quasipotential field U_A(x) using the Ordered
Line Integral Method (Dahiya & Cameron 2018). The state dimension D is taken
from sys::CoupledSDEs{IIP, D, I, P} and must match grid::CartesianGrid{D}.

# Keyword arguments
- band_radius::Int   = default_K(grid) -- K, the accepted-band radius in grid cells
                                          (default: `max(5, round(Int, sqrt(minimum(nbox))))`, capped at 32)
- near_source_layers::Int = 3      -- size of the analytic CARE seed box; 0 disables
- verbose::Bool      = false
- show_progress::Bool = true

A warning is emitted for D > 4; the algorithm still runs but per-axis grid
resolution typically becomes too coarse to be useful past D = 4.
"""
function quasipotential(sys::CoupledSDEs{IIP, D, I, P},
                        grid::CartesianGrid{D, T},
                        attractor::AbstractVector{<:Real};
                        band_radius::Int   = default_K(grid),
                        near_source_layers::Int = 3,
                        verbose::Bool      = false,
                        show_progress::Bool = true) where {IIP, D, I, P, T}
    length(attractor) == D || throw(DimensionMismatch(
        "attractor has length $(length(attractor)) but sys has D=$D"))
    x_A = SVector{D, T}(attractor)
    return _quasipotential_impl(sys, grid, x_A,
                                Val(band_radius), Val(near_source_layers),
                                verbose, show_progress)
end
```

Exports: `quasipotential`, `QuasiPotential`.

The public function accepts any `AbstractVector{<:Real}` for ergonomics (users typically pass `[x, y]`), converts it once to `SVector{D, Float64}`, and wraps the `Int` kwargs into `Val{band_radius}`, `Val{near_source_layers}` before calling the internal type-stable `_quasipotential_impl`. All inner functions specialize on those `Val` types and on the `SVector` element type.

## 4. Algorithm

OLIM is a Dijkstra-style one-pass sweep of the static Hamilton-Jacobi equation $H(x, \nabla U) = 0$ with $H(x, p) = p \cdot b(x) + \tfrac{1}{2} p \cdot Q(x) p$, where $Q(x) = a(x)/s$ is the trace-normalized diffusion tensor matching `fw_action` (`src/largedeviations/action.jl`), $a = \sigma \sigma^\top$, and $s = \operatorname{tr}(a(u_0))/D$ pinned at the attractor.

The quantity propagated along straight segments is the **geometric Lagrangian**

$$L_g(x, v) \;=\; \|v\|_Q \, \|b(x)\|_Q \;-\; \langle v, b(x)\rangle_Q, \qquad \|w\|_Q^2 := w^\top Q(x)^{-1} w,$$

which is the Maupertuis (arc-length) reduction of $\tfrac{1}{2}\|v - b\|_Q^2$ after minimization over time parametrization. $L_g \ge 0$ with equality on the deterministic flow, and $L_g$ is monotone in path length. These are the properties that make Dijkstra-style one-pass sweeping correct on a static H-J. The naive time-parametrized $\tfrac{1}{2}\|v-b\|_Q^2$ does **not** have this monotonicity and cannot be used.

Using trace-normalized $Q$ ensures that `quasipotential(sys, grid, x_A)` agrees with `fw_action(sys, ...)` minimized over paths reaching the same endpoint. The trace factor $s$ is obtained via `CriticalTransitions._trace_normalized_a(sys)` once at construction.

**Near-zero drift.** $L_g$ has a removable singularity at $b = 0$ (the integrand goes to zero, but $\partial L_g / \partial v$ has a $1/\|b\|_Q$ term). Inside a tolerance ball $\|b\|_Q \le \varepsilon_b$ (default $\varepsilon_b = 10^{-10}$) the integrand is replaced by its Taylor expansion $L_g \approx \tfrac{1}{2}\|v\|_Q^2 - \langle v, b\rangle_Q + O(\|b\|_Q^2)$ and the closed-form gradient by its analytic limit. This matters since the sweep starts at the attractor where $b = 0$.

### 4.1 Outer loop

Four node statuses: `UNKNOWN` (no tentative value), `CONSIDERED` (in the heap, tentative U), `FRONT` (accepted, has at least one non-accepted neighbor), `ACCEPTED` (accepted, fully interior). Status is stored as `Array{UInt8, D}`. The set of FRONT cells is mirrored as a `BitArray{D}` for O(1) membership tests; simplex update enumeration iterates only the FRONT cells inside the K-stencil, not the full accepted region.

```
seed_near_source!(state, Val(D), Val(K_seed))   # analytic CARE box, status = FRONT
push every CONSIDERED cell adjacent to the seed box onto heap

while heap not empty:
    pop (_, lin_c) -> c
    status[c] = FRONT
    front[c]  = true

    for x in circular_stencil(c, Val(K), bounds):   # compile-time circular K-ball
        status[x] == ACCEPTED && continue
        (newU, newback) = local_update(Val(D), x, state, L, Val(K))
        newU < U[x] || continue
        U[x] = newU
        back_pointer[x] = newback
        if status[x] == UNKNOWN
            status[x] = CONSIDERED
            handles[lin(x)] = push!(heap, (newU, lin(x)))
        else                                         # already CONSIDERED
            update!(heap, handles[lin(x)], (newU, lin(x)))
        end
    end

    # Prune c from FRONT if no neighbor remains non-accepted.
    if all(status[n] >= FRONT for n in neighbors(c)):
        status[c] = ACCEPTED
        front[c]  = false
    end
end

return QuasiPotential(U, back_pointer, source, grid)
```

Two paths into the heap: `push!` for `UNKNOWN -> CONSIDERED` transitions, `update!` (decrease-key) for already-`CONSIDERED` cells whose tentative value improved.

**Heap.** `DataStructures.MutableBinaryHeap{Tuple{T, Int}, FasterForward}` plus a `Vector{Int}` of size $N^D$ (not a `Dict`) for linear-cell-index to heap-handle mapping. The vector is O(N^D) memory but constant-time lookup, and N^D is already paid by `U` and `status`.

**Stencil.** Compile-time **circular** stencil generated by a `@generated function stencil_offsets(::Val{K}, ::Val{D})` returning an `NTuple{N, CartesianIndex{D}}` of offsets with $\|\delta\|_2 \le K$. The geometric Lagrangian only cares about the line direction, so directional coverage is what matters; the circular stencil cuts the offset count by ~$\pi/4$ in 2D and ~$\pi/6$ in 3D vs. the hyperbox.

### 4.2 Local update

For a Considered cell $x$:

1. Iterate the precomputed circular stencil around $x$. For each offset, check `status` and short-circuit if not `FRONT` or `ACCEPTED` with U finite.
2. **Vertex candidates** (all in-range FRONT or close ACCEPTED cells $y$): $\Phi_v(y) = U(y) + \int_0^1 L_g(y + \tau v,\, v)\, d\tau$, $v = x - y$.
3. **Simplex candidates** (D=2: front edges; D=3: front triangles, dispatched on `Val{D}`, no-op otherwise): closed-form root-find on the KKT conditions (§4.4, §4.5).
4. Return the overall minimum and a `BackRef{D}` encoding the winning simplex.

The loop is allocation-free: $L_g$ evaluations use stack-allocated `SVector` arithmetic, the candidate `Φ`/back-ref pair is kept in a `Tuple{T, BackRef{D}}` accumulator. CI asserts `@allocated local_update(...) == 0` for a representative call.

Per-$D$ specialization is achieved by method dispatch on `Val{D}`:

```julia
add_simplex_candidates!(best, x, accepted, state, L, ::Val{D}) where {D} = best
add_simplex_candidates!(best, x, accepted, state, L, ::Val{2}) = ...
add_simplex_candidates!(best, x, accepted, state, L, ::Val{3}) = ...
```

The default method is a no-op (vertex-only path); specialized methods exist only for $D = 2$ and $D = 3$. The compiler picks the right one at each call site; there is no runtime branch on $D$.

### 4.3 Vertex candidate

Straight segment from accepted $y$ to Considered $x$, $v = x - y$:

$$\Phi_v(y) \;=\; U(y) \;+\; \int_0^1 L_g\big(y + \tau v,\; v\big)\, d\tau,\qquad L_g(z, v) = \|v\|_{Q(z)} \, \|b(z)\|_{Q(z)} - \langle v, b(z)\rangle_{Q(z)}.$$

For state-dependent $Q$, the metric is re-evaluated at every quadrature node along the segment (see §5). The integration variable is $\tau$ to avoid clashing with the trace factor $s$ of §4.1.

### 4.4 Edge minimization (D = 2)

For an admissible pair $(y_0, y_1)$ of mutually adjacent FRONT vertices, parametrize the back-point as $y(\lambda) = (1-\lambda)\,y_0 + \lambda\, y_1$, $\lambda \in [0, 1]$:

$$\Phi(\lambda) \;=\; U_{\text{interp}}(\lambda) \;+\; \int_0^1 L_g\big(y(\lambda) + \tau (x - y(\lambda)),\; x - y(\lambda)\big)\, d\tau.$$

**`U` interpolation.** Hermite cubic in $\lambda$ using $U(y_0)$, $U(y_1)$, and one-sided gradient estimates $\partial_\lambda U(y_0)$, $\partial_\lambda U(y_1)$ obtained from finite differences against the next-adjacent FRONT cell along the edge direction (skipped to linear if no such cell exists; degrades quietly to 2nd order). This matches the 4th-order accuracy of Simpson on $L_g$.

**Solver: root-find, not minimize.** For Simpson quadrature in $\tau$ with fixed nodes $\{0, 1/2, 1\}$, $\partial \Phi / \partial \lambda$ has a closed form (chain rule through $L_g$ and `U_interp`). Solve $\partial \Phi / \partial \lambda = 0$ on $[0, 1]$ via a **hand-rolled ITP** (Interpolate-Truncate-Project, Oliveira-Takahashi 2020) in `quasipotential/itp.jl`, ~40 lines, allocation-free, type-parametric on the scalar type `T`:

```julia
function itp_root(f, a::T, b::T;
                  k1 = T(0.1), k2 = T(2), n0 = 1,
                  atol = eps(T)^(T(3)/T(4)),
                  maxiter = 64) where {T <: AbstractFloat}
    # returns (root::T, converged::Bool); converged=false iff no sign change.
end
```

ITP has provably better worst-case complexity than Brent/Dekker (minmax-optimal between bisection and superlinear). Bracket check `f(a)*f(b) > 0` short-circuits to `(a or b, false)`; the caller treats that as "no interior minimum" and falls back to the endpoints. No new dependency, fully inlined into the hot loop.

Workflow per edge candidate: evaluate $\Phi(0)$, $\Phi(1)$, run `itp_root(dΦdλ, 0, 1)`, and if it converges with $\lambda^* \in (0,1)$ also evaluate $\Phi(\lambda^*)$. Take the smallest of the (up to three) candidate values.

The closed-form $\partial L_g / \partial \lambda$ uses the near-zero-drift branch of §4 whenever $\|b\|_Q \le \varepsilon_b$ at any quadrature node along the segment.

### 4.5 Triangle minimization (D = 3)

For an admissible triple $(y_0, y_1, y_2)$ of mutually adjacent FRONT vertices, parametrize the back-point in barycentric coordinates:

$$y(\lambda) = (1 - \lambda_1 - \lambda_2) y_0 + \lambda_1 y_1 + \lambda_2 y_2,\qquad \lambda \in \Delta_2 := \{\lambda_1, \lambda_2 \ge 0,\ \lambda_1 + \lambda_2 \le 1\}.$$

**Solver: active-set KKT.** With Simpson on $L_g$ and bilinear $U$-interpolation on the triangle, $\nabla \Phi$ and $\nabla^2 \Phi$ have closed forms. Procedure (per candidate triangle):

1. Take a Newton step on the unconstrained gradient from the centroid.
2. If the iterate lies in $\Delta_2$ and $\|\nabla \Phi\|_\infty < \mathrm{tol}$, accept the interior minimum.
3. Otherwise activate the violated constraint(s) and reduce to:
   - one edge (1D root-find on $\partial \Phi / \partial t$ along the edge, same as §4.4),
   - or a vertex (no further minimization, just evaluate).
4. Take the smallest of the interior candidate (if found) and the three edge/vertex restrictions.

No outer-loop optimizer (no Nelder-Mead, no grid scan), no allocation. The active-set logic is `~80` lines and the closed-form gradient/Hessian are derived once in the docstring for `triangle_minimum`.

The three constituent edges are *always* evaluated (cheap, three Brent root-finds), so even a degenerate Newton step does not miss a boundary minimum.

### 4.6 Admissibility of simplices

A pair (or triple) is admissible iff:

- all vertices have `status == FRONT` (the explicit front set, not full ACCEPTED),
- they are mutually adjacent on the grid (Chebyshev distance 1 between every pair),
- they lie within the circular K-stencil of $x$,
- $x$ projects strictly outside the simplex on the "Unknown" side (i.e. the line from the closed-form minimizer to $x$ does not re-cross the front).

The admissibility check is per-D via `front_edges(x, Val(K))` / `front_triangles(x, Val(K))` iterators that read from the `BitArray{D}` front mask. With the front pruned (§4.1), the iterator yields O(K^(D-1)) candidates per call, not O(K^D).

### 4.7 Back-pointer convention

The `BackRef{D}` per cell stores `(v0, v1, s)`:
- one-point (vertex) update: `v1 == v0`, `s = NaN32`. Predecessor is exactly `v0`.
- two-point (edge) update: `v0`, `v1` are the adjacent FRONT vertices, `s ∈ [0,1]` (as `Float32`) is the interior parameter. Predecessor is `(1 - s)*v0 + s*v1`.
- three-point (triangle, D=3) update: encoded as two consecutive edge entries (`v0=y0`, `v1=barycenter-on-edge-y1-y2`), or by storing a small auxiliary `Dict{CartesianIndex{3}, NTuple{3, Float32}}` for the few cells actually hit by an interior-triangle minimizer. (Triangle interior minima are rare in practice; this is a minor extension over the dense `BackRef` array.)

Instanton reconstruction walks the chain: at each step, decode `BackRef`, take the sub-cell continuous predecessor, and find the cell containing it; repeat until reaching the source cell. This yields a discrete-but-sub-cell instanton without needing `minimize_geometric_action` for refinement (still available as an optional polish step).

## 5. Quadrature

Simpson's rule along the unit-parameter segment, inlined at every call site:

```julia
@inline line_integral(L::GeometricLagrangian{D, T}, y, v) where {D, T} =
    (L(y, v) + 4 * L(y + T(0.5) * v, v) + L(y + v, v)) / 6
```

Three evaluations of $L_g$ per candidate segment. Each evaluation re-computes $Q(z)^{-1}$ at $z$ (a no-op for additive noise since $Q_{\mathrm{inv}}$ is then a stored `SMatrix`; a fresh `inv(SMatrix)` for multiplicative noise). Quadrature order matches the §4.4 Hermite-cubic `U`-interpolation, so the full vertex/edge update is 4th-order accurate in the step length when both $b$ and $Q$ are smooth.

## 6. Geometric Lagrangian wrapper

Internal type, not exported; binds the drift, the trace-normalized inverse metric, and a near-zero-drift threshold:

```julia
struct GeometricLagrangian{D, T, B, A}
    b::B          # SVector{D, T} -> SVector{D, T}
    Q_inv::A      # SMatrix{D,D,T,D*D} for additive, x -> SMatrix for multiplicative
    eps_b::T      # near-zero-drift threshold (Q-norm), default T(1e-10)
end

@inline function (L::GeometricLagrangian{D, T})(x::SVector{D, T},
                                                 v::SVector{D, T}) where {D, T}
    bx   = L.b(x)
    Qinv = _Q_inv_at(L.Q_inv, x)
    bnQ2 = dot(bx, Qinv, bx)        # |b|_Q^2
    if bnQ2 < L.eps_b * L.eps_b
        # Taylor branch (removable singularity), see §4 near-zero drift.
        return T(0.5) * dot(v, Qinv, v) - dot(v, Qinv, bx)
    end
    vnQ2 = dot(v,  Qinv, v)
    return sqrt(vnQ2) * sqrt(bnQ2) - dot(v, Qinv, bx)
end

@inline _Q_inv_at(M::SMatrix, _)              = M
@inline _Q_inv_at(f::F, x) where {F<:Function} = f(x)
```

Construction:

* drift via `dynamic_rule(sys)`. For an out-of-place `CoupledSDEs` the function is used directly; for an in-place one it is wrapped in a closure that fills a stack-allocated `MVector{D, T}` and returns an `SVector{D, T}`. The `IIP` flag is consumed at construction time and does not appear on the struct type.
* $Q$ via `CriticalTransitions._trace_normalized_a(sys)` (`src/largedeviations/utils.jl:144`). For additive noise (`Base.Returns`) the constant matrix is inverted once into `SMatrix{D, D, T, D*D}` and stored; for state-dependent noise it is wrapped as `x -> inv(SMatrix{D, D, T}(Q(x)))`.

The struct is internal-only (no export, no compat guarantee). It exists to bundle the four pieces (drift, metric, magnitude tolerance, dispatch type) that every `L_g` evaluation needs, and to keep the signature of `line_integral`, `vertex_candidate`, `edge_minimum`, `triangle_minimum` to a single Lagrangian argument.

## 7. Near-source analytic seeding

Inside a box of side $2 K_{\text{seed}} + 1$ centered at the attractor cell, set $U(x) \approx (x - x_A)^\top P (x - x_A)$ where $P$ solves the continuous-time algebraic Riccati equation

$$A^\top P + P A + 2\, P\, Q\, P = 0,\qquad A = \nabla b(x_A),\ Q = a(x_A)/s.$$

(Same trace-normalized $Q$ as §4.1.) Solver: `MatrixEquations.arec` for $D \geq 3$, closed-form for $D = 2$. The drift Jacobian $A$ is obtained via `ForwardDiff.jacobian` on the OOP drift wrapper built in §6.

Cells inside the box are marked `ACCEPTED` with their analytic $U$ values. The first sweep iteration picks up from the box edge.

Guard: if `eigvals(A)` has any non-negative real part, $x_A$ is not a stable fixed point of $b$; skip analytic seeding (emit `@info` when `verbose`), start the sweep from the source cell with $U = 0$ and `near_source_layers = 0` semantics.

`near_source_layers = 0` is supported and disables analytic seeding.

## 8. Compile-time dispatch summary

Every dimension-dependent code path is selected by Julia method dispatch:

| Decision | Dispatched on |
|---|---|
| state dimension | `Val{D}` (specialized methods for `Val{2}`, `Val{3}`, generic fallback) |
| band radius (stencil size) | `Val{K}` |
| near-source box size | `Val{K_{\text{seed}}}` |
| in-place vs out-of-place drift | `IIP` parameter inherited from `CoupledSDEs` |
| constant vs state-dependent noise | type of `GeometricLagrangian.Q_inv` (`SMatrix` vs callable) |

There is no `if D == ...` or `if band_radius == ...` branch at runtime inside any hot loop.

## 9. File layout

```
src/largedeviations/quasipotential/
├── quasipotential.jl    public API, QuasiPotential, BackRef, _quasipotential_impl, default_K
├── state.jl             solver state container (U, back_pointer, status, front, heap, handles)
├── stencil.jl           @generated circular K-stencil, neighbor iterators
├── front.jl             FRONT bit-mask maintenance and front_edges / front_triangles iterators
├── sweep.jl             Dijkstra sweep loop, front pruning
├── update.jl            local_update dispatcher, vertex_candidate
├── simplex.jl           edge_minimum (D=2, ITP root-find on dΦ/dλ), triangle_minimum (D=3, active-set KKT)
├── itp.jl               hand-rolled ITP scalar bracketed root-find (~40 lines)
├── quadrature.jl        line_integral (Simpson), hermite_U_interp
├── lagrangian.jl        GeometricLagrangian struct and constructors (Taylor branch for |b| ~ 0)
└── source_seed.jl       CARE solve (MatrixEquations.arec / closed-form 2D), seed_near_source!
```

`src/CriticalTransitions.jl` adds `include("largedeviations/quasipotential/quasipotential.jl")` and exports `quasipotential`, `QuasiPotential`, `BackRef`. `quasipotential/quasipotential.jl` includes its siblings in dependency order.

## 10. Dependencies

Verified against `Project.toml` on `main` (2026-05-25):

| Package | Status | Compat | Reason |
|---|---|---|---|
| `DataStructures` | **add** | ≥ 0.18 | `MutableBinaryHeap{Tuple{T, Int}, FasterForward}` for the Dijkstra priority queue with decrease-key |
| `MatrixEquations` | **add** | ≥ 2.0 | `arec` for the CARE-based near-source seed |
| `ForwardDiff` | present | (existing) | Jacobian of drift at the attractor |
| `ProgressMeter` | present | (existing) | `Progress` / `next!` for `show_progress=true` |
| `StaticArrays`, `LinearAlgebra` | present | (existing) | `SVector`/`SMatrix`, `dot`, `inv` |
| `DynamicalSystemsBase` | present | (existing) | `CoupledSDEs`, `dynamic_rule`, `dimension`, `diffusion_function` |

Explicitly not used (hand-rolled instead, with rationale):

* `Roots.jl` / `BracketingNonlinearSolve.jl` for the 1D edge root-find. A ~40-line hand-rolled ITP (§4.4) is allocation-free, type-parametric on `T`, and avoids per-call problem/solution boilerplate in the SciML option. The maintenance cost is negligible since ITP is short and the algorithm is closed.
* `Optim.jl` / `OptimizationBase.jl` for the 3D triangle update. Active-set KKT with closed-form gradient/Hessian (§4.5) is `~80` lines, allocation-free, and avoids the overhead of a generic constrained optimizer.
* `Interpolations.jl` / `FastInterpolations.jl` for the Hermite cubic along edges. Closed-form (§4.4) is 5 lines and fully unrolls.
* `Stencils.jl` / `OffsetArrays.jl` for the K-stencil. A `@generated` function returning `NTuple{N, CartesianIndex{D}}` is cleaner and compile-time known.
* `HCubature.jl` / `QuadGK.jl` for the line integral. Fixed 3-node Simpson (§5) is the right scope.
* `NonlinearSolveBase.jl` for the root-find. Overkill for a 1D scalar problem; `Roots.A42` is the right tool.

## 11. Testing

`test/largedeviations/quasipotential.jl`:

1. **Geometric Lagrangian unit** - $L_g$ matches $\|v\|_Q \|b\|_Q - \langle v, b\rangle_Q$ to machine precision for additive and multiplicative noise; Taylor branch agrees with the full formula to `1e-12` at $\|b\|_Q = 10\,\varepsilon_b$.
2. **Vertex update unit** - 2D quadratic well, analytic $U$; verify single-cell `vertex_candidate` to expected order.
3. **Edge root-find unit (D=2)** - synthetic problem with analytic optimal $\lambda^*$; the ITP root-find on $\partial \Phi / \partial \lambda$ recovers $\lambda^*$ to `1e-10` in `≤ 12` iterations, and returns `(_, false)` cleanly on a monotone test case.
4. **Triangle KKT unit (D=3)** - synthetic with known interior minimum; verify active-set selects the correct constraint face for boundary cases.
5. **Front pruning** - after sweep on a $50^2$ grid, the FRONT bit-mask contains only cells with at least one non-accepted neighbor; ACCEPTED interior count is consistent with the total minus FRONT.
6. **Allocation** - `@allocated local_update(...)` is 0 for a representative call in D=2 and D=3, additive and multiplicative noise.
7. **End-to-end gradient (2D double well)** - $U = 2V$ on a $100^2$ grid; relative error tolerance 1%.
8. **End-to-end Maier-Stein (2D non-gradient)** - $U$ at the saddle vs the literature value (≈ 0.5); also vs `minimize_geometric_action`. Tolerance 2%. The OLIM value must agree with `fw_action(sys, instanton, time)` for the gMAM instanton (both use the trace-normalized $Q$ convention).
9. **3D linear Ornstein-Uhlenbeck smoke** - analytic quadratic $U$; coarse $30^3$ grid; mid-grid relative error < 5%.
10. **D=5 warning** - tiny $15^5$ grid; `@test_logs (:warn, r"quasipotential in D=5") quasipotential(...)`.
11. **Back-pointer walk** - from the saddle cell, follow the `BackRef` chain (decoding `s` for two-point updates); verify termination at the source, monotone-decreasing $U$, and that the sub-cell instanton has lower discrete action than the nearest-vertex one.
12. **Multiplicative-noise (2D)** - Maier-Stein drift with state-dependent diagonal $\sigma(x) = \mathrm{diag}(1, 1 + 0.3 x_1^2)$; solver runs to completion, agrees with `fw_action(sys, instanton, time)` on the gMAM instanton to 5%, `_trace_normalized_a(sys)` is threaded through to $L_g$ correctly.

Tests follow the existing `test/largedeviations/` style: `@testset` blocks, seeded RNG, `CoupledSDEs` constructors with out-of-place drifts returning `SVector`.

## 12. Example

`examples/quasipotential_maierstein.jl`:

- Define the Maier-Stein system as a `CoupledSDEs` (matches `examples/gMAM_Maierstein.jl`).
- Build a $200^2$ `CartesianGrid` covering $[-1.5, 1.5]^2$.
- Run `quasipotential(sys, grid, SVector(-1.0, 0.0))` from the left attractor.
- Visualize $U$ with `CairoMakie.contourf`, overlay the gMAM instanton, mark attractor and saddle.
- Print the OLIM action at the saddle cell and the gMAM action for comparison.

## 13. Documentation

- `docs/src/man/largedeviations.md` (or the analogous existing page) gains a section on `quasipotential` (briefly noting OLIM as the algorithm) with a one-screen usage block and a pointer to the Maier-Stein example.
- `docs/src/refs.bib` gains the Dahiya-Cameron 2018 entry and Cameron's earlier 2012 *Physica D* paper.
- `CHANGELOG.md` gets an entry under the next release.

## 14. Risks and known limitations

- **Memory** scales as $N^D$ (U, status, front bitmask, BackRef array, heap-handle vector). The `D > 4` warning is informational; we do not refuse to run.
- **Heap operations** are $O(\log N^D)$ per pop/push; the dominant cost in practice is the simplex update inside the K-stencil, not the heap.
- **State-dependent noise** is benchmarked on the 2D smoke test (test 12) but not on a published reference solution; expect to revisit accuracy when one becomes available.
- **Caustics** are not detected. The sweep returns the viscosity solution (correct minimum), but users wanting to visualize Maxwell sets will need to layer that on themselves in a follow-up.
- **Trace-normalization pin point.** $s = \operatorname{tr}(a(u_0))/D$ is pinned at the attractor cell. For state-dependent noise this matches `fw_action`, and means quasipotential values are comparable across calls only when `x_A` is the same.
- **K too small.** `default_K(grid)` is a heuristic. Users with very anisotropic problems (slow drift in one direction) may need to tune K manually; we document this in the example.

## 15. References

1. S. Dahiya and M. K. Cameron. *Ordered Line Integral Methods for Computing the Quasi-Potential.* Journal of Scientific Computing 75 (2018). arXiv:1808.00562. Cameron's accompanying code: https://math.umd.edu/~mariakc/OLIM.html
2. M. K. Cameron. *Finding the quasipotential for nongradient SDEs.* Physica D 241 (2012).
3. J. A. Sethian and A. Vladimirsky. *Ordered upwind methods for static Hamilton-Jacobi equations.* SIAM Journal on Numerical Analysis 41 (2003).
4. T. Grafke, T. Schäfer, and E. Vanden-Eijnden. *Long term effects of small random perturbations on dynamical systems: theoretical and computational tools.* In *Recent Progress and Modern Challenges in Applied Mathematics* (2019); arXiv:1604.03818.
5. N. Yang, S. F. Potter, and M. K. Cameron. *Computing the quasipotential for nongradient SDEs in 3D.* Journal of Computational Physics 379 (2019). Source of the active-set KKT update for the D=3 triangle minimization (§4.5).
6. I. F. D. Oliveira and R. H. C. Takahashi. *An Enhancement of the Bisection Method Average Performance Preserving Minmax Optimality.* ACM Transactions on Mathematical Software 47:1 (2020). ITP algorithm used in §4.4 (hand-rolled in `itp.jl`).
