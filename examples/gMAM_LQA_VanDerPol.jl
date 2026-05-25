# # gMAM with Local Quadratic Approximation: escape from a limit cycle
#
# We compute the minimum-action escape path from a stable limit cycle to a saddle point in a planar SDE, following Section 6.2 of [lin_quasi_2018](@citet). The target value is `V ≈ 0.1567` from [cameron_finding_2012](@citet) for the forced Van der Pol oscillator.
#
# The SDE `dX_t = b(X_t)\,dt + \sqrt{ε}\,σ(X_t)\,dW_t` has a stable limit cycle `Γ` and an unstable fixed point `x_s` inside it. As `ε → 0` the mean escape time from `Γ` over `x_s` scales as `\exp(V(x_s)/ε)`, where `V` is the Freidlin-Wentzell quasi-potential `V(x) = \inf_{T,φ:Γ→x} S[φ]` with action `S[φ] = \tfrac12 \int (\dot φ - b)^\mathsf{T} a^{-1} (\dot φ - b)\,dt` and `a = σ σ^\mathsf{T}`.
#
# Computing `V` directly is hard because `S` must be minimized over both the path and its duration. Near `Γ` however `V` is smooth and quadratic in the transverse direction, so we can pre-compute its local shape and use it to anchor a global path optimization. This is gMAM-LQA: build a moving frame along `Γ`, solve a periodic Riccati equation for the transverse Hessian `G(τ)` of `V`, then run gMAM outside the tube with the launch point on `Γ` optimized by projected gradient on the augmented action.

using CriticalTransitions
using PeriodicOrbits
using Attractors
using CairoMakie
using LinearAlgebra: norm

# The forced Van der Pol oscillator,
# ```math
# \dot{x} = x - \tfrac{1}{3} x^3 + y - \tfrac{1}{9} y^3, \qquad \dot{y} = x + 0.9 ,
# ```
# has one stable limit cycle and one unstable fixed point at `x_s ≈ (-0.9, 0.6942)`. The forcing `+0.9` breaks the standard Van der Pol symmetry.

function vdp(u, p, t)
    x, y = u
    return SA[x - x^3 / 3 + y - y^3 / 9, x + 0.9]
end

sys_det = CoupledODEs(vdp, [2.0, 0.0])

# Find the limit cycle with `AttractorsViaRecurrences`. We only need one point on `Γ`; the recurrent set returned by the mapper is too coarse to discretize the periodic Riccati equation, so we use it just as a seed.

grid = (range(-3, 3; length = 81), range(-4, 4; length = 81))
mapper = AttractorsViaRecurrences(sys_det, grid;
    sparse = true, consecutive_recurrences = 400,
)
sampler, _ = statespace_sampler(HRectangle([-2.5, -3.5], [2.5, 3.5]), 42)
fractions = basins_fractions(mapper, sampler; N = 30)
attractors = extract_attractors(mapper)

# `PeriodicOrbits.PeriodicOrbit(ds, u0, T)` integrates the deterministic flow for exactly one period at the resolution we want. The period `T ≈ 5.7966` is the value from Section 6.2; in general it can be obtained from `OptimizedShooting` or a Poincaré-section analysis.

u_on_Γ = collect(first(attractors[1]))
T_period = 5.7966
Nτ = 200
po = PeriodicOrbit(sys_det, u_on_Γ, T_period; Δt = T_period / Nτ)

# The noise strength is set to one; the trace-normalized action used internally is invariant under rescalings of `σ`. `LimitCycleFrame(po, sys)` builds the moving frame: an orthonormal `d × d` basis `E(τ)` at each orbit sample, with first column tangent to `Γ` and remaining columns spanning the transverse subspace (the latter stored separately as `Ẽ`). It also computes the projected Jacobian `M̃ = Ẽ^\mathsf{T} (J - \dot E E^\mathsf{T}) Ẽ` and the trace-normalized projected noise covariance `Ã`. The frame is transported by the variational equation `\dot Φ = J(γ(τ))\,Φ` and sign-continued to avoid Floquet flips.

σ = 1.0
sys = CoupledSDEs(vdp, u_on_Γ; noise_strength = σ)
lc = LimitCycleFrame(po, sys)

# In a tube of radius `h` around `Γ`, writing `x = γ(τ) + Ẽ(τ)\,z`, the quasi-potential admits the expansion `V(x) = \tfrac12 z^\mathsf{T} G(τ) z + O(\|z\|^3)`, where `G(τ)` is the unique positive-definite `T`-periodic solution of
# ```math
# \dot G = -\tilde M^\mathsf{T} G - G \tilde M - G \tilde A G .
# ```
# `local_quasipotential(lc)` solves this fixed-point equation by repeated integration over one period until the period map converges. In two dimensions `G(τ)` is a scalar.

G = local_quasipotential(lc)

# Outside the tube we minimize the geometric action `\hat S[φ]` over arclength-parametrized paths with prescribed endpoints. The first endpoint, however, is *not* fixed: only the tube radius `h` is constrained, and the launch parameters `(τ_\text{start}, \hat z)` are optimized jointly with the path on the augmented action
# ```math
# \hat S_\text{tot} = \hat S[φ] + \tfrac{1}{2} h^2 \hat z^\mathsf{T} G(τ_\text{start}) \hat z .
# ```
# The algorithm alternates inner gMAM sweeps with frozen endpoints and outer projected-gradient steps on `(τ_\text{start}, \hat z)`. The converged value of `\hat S_\text{tot}` equals the quasi-potential at the target, up to an `O(h^3)` discretization error.

x_saddle = SA[-0.9, 0.6942]
mp = minimize_geometric_action(
    sys, lc, x_saddle;
    G = G,
    tube_radius = 0.05,
    npoints = 160,
    N_inner = 300,
    maxiters = 200,
    reltol = 1.0e-7,
    stepsize = 0.005,
)

# Below we plot the optimal escape path together with `G(τ)`, the transverse Hessian of `V` along `Γ`. The launch point on `Γ` is marked with a gold star, the saddle with a white circle.

γ_mat = reduce(hcat, [collect(p) for p in lc.γ])
G_scalar = vec(G[1, 1, :])
map_pts = reduce(hcat, [collect(p) for p in mp.path])

x_lo, x_hi = -2.5, 1.5
y_lo, y_hi = -0.5, 4.2
odeSol(x, y) = Point2f(vdp(SA[x, y], nothing, 0.0)...)
launch_idx = argmin([norm(γ_mat[:, k] .- map_pts[:, 1]) for k in 1:Nτ])

set_theme!(theme_light(); fontsize = 13, fonts = (regular = "TeX Gyre Heros", italic = "TeX Gyre Heros Italic"))

fig = Figure(; size = (1000, 460))

ax1 = Axis(fig[1, 1];
    xlabel = "x", ylabel = "y",
    title = "Optimal escape path  Γ → saddle",
    aspect = DataAspect(),
)
streamplot!(ax1, odeSol, (x_lo, x_hi), (y_lo, y_hi);
    gridsize = (22, 22), arrow_size = 7, stepsize = 0.01,
    colormap = [(:gray60, 0.5), (:gray60, 0.5)],
)
lines!(ax1, γ_mat[1, :], γ_mat[2, :];
    color = :crimson, linewidth = 3, label = "limit cycle Γ",
)
lines!(ax1, map_pts[1, :], map_pts[2, :];
    color = :darkorange, linewidth = 3, label = "MAP",
)
scatter!(ax1, [γ_mat[1, launch_idx]], [γ_mat[2, launch_idx]];
    marker = :star5, color = :gold, strokecolor = :black, strokewidth = 1.2,
    markersize = 18, label = "launch point",
)
scatter!(ax1, [x_saddle[1]], [x_saddle[2]];
    color = :white, strokecolor = :black, strokewidth = 1.4, markersize = 14,
    label = "saddle",
)
axislegend(ax1; position = :rb, framevisible = true, padding = 6, rowgap = 1)
limits!(ax1, x_lo, x_hi, y_lo, y_hi)

ax2 = Axis(fig[1, 2];
    xlabel = "τ / T",
    ylabel = "transverse Hessian  G(τ)",
    title = "Quasi-potential coefficient along Γ",
)
τ_norm = range(0, 1; length = Nτ)
band!(ax2, τ_norm, zero(τ_norm), G_scalar; color = (:indigo, 0.15))
lines!(ax2, τ_norm, G_scalar; color = :indigo, linewidth = 2.5)
vlines!(ax2, [launch_idx / Nτ]; color = :gold, linewidth = 2, linestyle = :dash,
    label = "τ_start / T",
)
axislegend(ax2; position = :rt, framevisible = false)
xlims!(ax2, 0, 1)
ylims!(ax2, 0, nothing)
fig

# The path leaves `Γ` from the side closest to the saddle and reaches `x_s` along the slow manifold. `G(τ)` is strongly modulated along the cycle, varying by more than an order of magnitude between its peak and its minimum on the slow branch.
#
# It is tempting to expect the launch to coincide with the minimum of `G(τ)`, since `\tfrac12 h^2 G(τ_\text{start})` is the LQA cost of crossing the tube boundary. It does not, and the reason is quantitative: at `h = 0.05` the boundary term is `\tfrac12 h^2 G \sim 10^{-3}` while the path action is `\sim 10^{-1}`, so the boundary contribution is sub-percent of the total. The minimization is dominated by `\hat S[φ]` and the launch ends up where `Γ` is closest to `x_s` in the Onsager-Machlup geodesic sense. The relative weight of the two terms is set by `h`: at smaller `h` the boundary term is even more negligible, while at larger `h` it starts to bias the launch toward small `G(τ)`, at the cost of degrading the LQA itself.
#
# A second view of the same information, with `Γ` colored by `G(τ)`:

fig2 = Figure(; size = (620, 480))
ax = Axis(fig2[1, 1];
    xlabel = "x", ylabel = "y",
    title = "Local quasi-potential stiffness G(τ) on Γ",
    aspect = DataAspect(),
)
streamplot!(ax, odeSol, (x_lo, x_hi), (y_lo, y_hi);
    gridsize = (22, 22), arrow_size = 7, stepsize = 0.01,
    colormap = [(:gray70, 0.5), (:gray70, 0.5)],
)
sc = scatter!(ax, γ_mat[1, :], γ_mat[2, :];
    color = G_scalar, colormap = :viridis, markersize = 8,
)
Colorbar(fig2[1, 2], sc; label = "G(τ)")
lines!(ax, map_pts[1, :], map_pts[2, :]; color = :darkorange, linewidth = 3)
scatter!(ax, [γ_mat[1, launch_idx]], [γ_mat[2, launch_idx]];
    marker = :star5, color = :gold, strokecolor = :black, strokewidth = 1.2,
    markersize = 18,
)
scatter!(ax, [x_saddle[1]], [x_saddle[2]];
    color = :white, strokecolor = :black, strokewidth = 1.4, markersize = 14,
)
limits!(ax, x_lo, x_hi, y_lo, y_hi)
fig2

# The computed action is

mp.action

# against the [cameron_finding_2012](@citet) reference `0.1567`, with relative error

abs(mp.action - 0.1567) / 0.1567

# i.e. around half a percent, inside both the path discretization error and the `O(h)` correction expected from the LQA. The action value is the leading exponent of the Kramers escape time `\mathbb{E}\,τ \sim \exp(V/ε)`, so a small change in `V` is a multiplicative shift in the predicted escape time.
