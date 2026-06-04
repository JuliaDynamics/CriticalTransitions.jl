# Developer notes & internals

This page documents implementation internals that help contributors and explain performance characteristics, but are **not** part of the user-facing API contract. User-facing usage and gotchas live on the method manual pages, in particular [Large deviation theory](@ref).

## sgMAM: compile-time dispatch on the diffusion tensor

Compile-time dispatch on the shape of ``a(x)`` is driven by the type of `sys.a` (a `Base.Returns` wrapper marks constant-in-`x` diffusion) together with the type of `a(x_ref)` (a `LinearAlgebra.Diagonal` marks diagonal coupling). Classification happens once when the per-path cache is built; the resulting cache type then selects the diagonal or coupled inner loops in `update_p!`, `update_x!`, and `geometric_gradient_step!`.

Rank-deficient ``a(x)`` is rejected at cache build: `a` is probed at a reference state `x_ref` (the first path point) and at the ``2D`` neighbors `x_ref ± h·eₗ`.

## Multiple shooting: BVP internals

The boundary states `y_0, y_{N_seg}` are parameterized by the unstable / stable eigenvectors of the Hamiltonian Jacobian `M = [J A; 0 -Jᵀ]` at each fixed-point endpoint. The BVP is segmented into `nshoots` shooting segments stitched by Newton-iterated continuity, solved by `NonlinearSolveFirstOrder.NewtonRaphson` with a sparse `ForwardDiff` Jacobian (block-bidiagonal sparsity contracted by `SparseMatrixColorings.GreedyColoringAlgorithm`) and `LinearSolve.UMFPACKFactorization()`.

There is no runtime guard for interior fixed-point crossings: cheap detection requires integrating, which the residual function does each Newton iteration anyway. The user-facing consequence (split a through-saddle transition into uphill and downhill legs) is documented on the [Large deviation theory](@ref) page.