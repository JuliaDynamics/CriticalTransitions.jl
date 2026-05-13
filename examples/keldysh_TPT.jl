# # [Quasipotential and TPT for a Keldysh oscillator](@id keldysh_TPT_example)

# This example studies a stochastic oscillator with nonlinear damping that arises
# as the classical limit of the Keldysh field theory for a driven-dissipative
# quantum oscillator. The deterministic dynamics has a stable limit cycle around
# the origin and two stable fixed points at large amplitude. We use the layered
# API of `CriticalTransitions.jl` to:
#
# 1. find the deterministic attractors and basins (via `Attractors.jl`),
# 2. compute the small-noise quasipotential $\Phi \simeq -D\log\rho_{ss}$ from
#    the stationary Fokker–Planck solution on a Cartesian grid,
# 3. apply Transition Path Theory: forward committor, reactive density and
#    current, and the mean first-passage time.
#
# ## Noise-strength regimes
#
# The control parameter is $\Delta\Phi/D$, the barrier-to-noise ratio.
#
# - **Large-deviation regime ($\Delta\Phi/D \gg 1$).** Freidlin–Wentzell:
#   $\rho_{ss}\sim e^{-\Phi/D}$ concentrates on attractors, $\tau\sim e^{\Delta\Phi/D}$.
#   The transition geometry is $D$-independent.
# - **TPT regime (any $D>0$).** The committor's saddle layer has width $\sim\sqrt{D}$;
#   reactive currents spread over a tube. Natural for the grid-based solve below.
#
# **Practical limits for this PDE method.**
#
# | Constraint | Floor | Notes |
# |---|---|---|
# | spectral gap above `eps(T)` | $\Delta\Phi/D < \log(1/\varepsilon_T)$ | `Float64`: 36, `Double64`: 71, 256-bit `BigFloat`: 177 |
# | $\rho_{ss}$ not underflowing | $\Phi_{\max}/D < \log(1/\varepsilon_T)$ | same per-type bound |
# | saddle layer resolved by grid | $\sqrt{D} \gtrsim \Delta x$ | refine the grid for smaller $D$ |
#
# Grid discretisation smooths the true gap, so in practice `Double64` works
# well past its analytical floor (down to $D \sim 10^{-3}$ at the grid below);
# `BigFloat` is exact at any $D$ in the working range, ~6× slower. Below
# $D \approx 10^{-2}$ the geometry is already $D$-independent, so instanton
# solvers (`string_method`, gMAM) give the same answer far cheaper.
#
# ## Model
#
# The SDE comes from the Keldysh action
# ```math
# S = \int\mathrm{d}t\,\bigl[\dot{x}_{\mathrm{cl}}\dot{x}_{\mathrm{q}}
#     - \dot{x}_{\mathrm{cl}}x_{\mathrm{q}}\mathcal{D}(x_{\mathrm{cl}})
#     - x_{\mathrm{q}}V'(x_{\mathrm{cl}}) + iDx_{\mathrm{q}}^2\bigr],
# ```
# with non-linear damping and an anharmonic potential
# ```math
# \mathcal{D}(x) = \gamma - \mu(1-x^2),\qquad
# V'(x) = \omega_0^2\,x + \alpha_3\,x^3 + \alpha_5\,x^5 .
# ```
# Eliminating the quantum field $x_{\mathrm{q}}$ yields a second-order Langevin
# equation on the classical phase space $(x, p) \equiv (x_{\mathrm{cl}},
# \dot{x}_{\mathrm{cl}})$:
# ```math
# \dot{x} = p,\qquad
# \dot{p} = -p\,\mathcal{D}(x) - V'(x) + \sqrt{2D}\,\xi(t) ,
# ```
# with white noise on the momentum channel only — this matches the $iDx_{\mathrm{q}}^2$
# term in the action.

using CriticalTransitions, StaticArrays, ProgressMeter
import LinearSolve as LS
using Sparspak  # enables LS.SparspakFactorization() on non-Float64 sparse matrices
using Attractors
using LinearAlgebra: norm
using Statistics: quantile
using CairoMakie

const ω₀ = 1.0
const γ = 0.03
const μ = 0.1
const α₃ = -0.25
const α₅ = 0.01
const D = 0.001 # noise strength D in the action; momentum noise = √(2D)

𝒟(x) = γ - μ * (1 - x^2)
V′(x) = ω₀^2 * x + α₃ * x^3 + α₅ * x^5

drift_kel(u, p_, t) = SA[u[2], -u[2] * 𝒟(u[1]) - V′(u[1])]
diff_kel(u, p_, t) = SA[0.0, sqrt(2D)]

sys = CoupledSDEs(drift_kel, [0.5, 0.0]; g = diff_kel, noise_strength = 1.0)

# `covariance_matrix(sys)` is `Diagonal([0, 2D])` — rank-1 diffusion. The
# `DiffusionGenerator` discretisation handles this correctly: along the
# noise-free $x$ axis it falls back to a first-order upwind stencil.

# ## Deterministic phase portrait
#
# The fixed points come from $V'(x)=0$. With our parameters
# ```math
# V'(x) = x\bigl(1 - 0.25\,x^2 + 0.01\,x^4\bigr) ,
# ```
# the non-trivial roots satisfy $x^2 = 5$ or $x^2 = 20$. Linearising shows
# $(0, 0)$ is an unstable spiral, $(\pm\sqrt{5}, 0)$ are saddles, and
# $(\pm\sqrt{20}, 0)$ are stable spirals. Near the origin the Van-der-Pol
# pumping ($\mu > \gamma$) creates a stable limit cycle.

x_saddle = sqrt(5)
x_fp = sqrt(20)
@show x_saddle, x_fp

# Visualise the deterministic flow:

fig = Figure(; size = (760, 420))
ax = Axis(
    fig[1, 1]; xlabel = "x", ylabel = "p", aspect = DataAspect(),
    title = "Deterministic flow"
)
streamplot!(
    ax,
    u -> Point2f(u[2], -u[2] * 𝒟(u[1]) - V′(u[1])),
    -5.5 .. 5.5, -3.0 .. 3.0;
    linewidth = 0.7,
    colormap = [:gray40, :gray40],
    arrow_size = 7,
    gridsize = (32, 16),
)
scatter!(ax, [-x_fp, x_fp], [0.0, 0.0]; color = :red, markersize = 14, label = "stable FP")
scatter!(
    ax, [-x_saddle, x_saddle], [0.0, 0.0]; color = :orange,
    marker = :diamond, markersize = 12, label = "saddle"
)
axislegend(ax; position = :rb)
fig

# ## Attractors and basins
#
# We hand the deterministic system to `Attractors.jl`'s recurrence-based
# mapper. The same grid will reappear as the discretisation of the diffusion
# generator below.

ds = CoupledODEs(drift_kel, [0.5, 0.0])
xg = range(-5.5, 5.5; length = 251)
pg = range(-3.0, 3.0; length = 121)
mapper = AttractorsViaRecurrences(
    ds, (xg, pg);
    consecutive_recurrences = 200, Dt = 0.1, sparse = true
)
basins, attractors = basins_of_attraction(mapper, (xg, pg); show_progress = false)

# Identify the limit cycle as the attractor with the largest support, and the
# two stable fixed points by their location:

lc_key = argmax(k -> length(attractors[k]), keys(attractors))
lc_pts = [SA[u[1], u[2]] for u in attractors[lc_key]]

right_fp_key = argmin(
    k -> norm(
        [
            attractors[k][1][1] - x_fp,
            attractors[k][1][2],
        ]
    ),
    [k for k in keys(attractors) if k != lc_key]
)
left_fp_key = argmin(
    k -> norm(
        [
            attractors[k][1][1] + x_fp,
            attractors[k][1][2],
        ]
    ),
    [k for k in keys(attractors) if k != lc_key]
)

fig = Figure(; size = (760, 420))
ax = Axis(
    fig[1, 1]; xlabel = "x", ylabel = "p", aspect = DataAspect(),
    title = "Basins of attraction"
)
heatmap!(
    ax, xg, pg, basins;
    colormap = Makie.Categorical(:Set2)
)
scatter!(
    ax, [p[1] for p in lc_pts], [p[2] for p in lc_pts];
    color = :black, markersize = 3, label = "limit cycle"
)
scatter!(ax, [-x_fp, x_fp], [0.0, 0.0]; color = :red, markersize = 12)
scatter!(
    ax, [-x_saddle, x_saddle], [0.0, 0.0]; color = :orange,
    marker = :diamond, markersize = 10
)
axislegend(ax; position = :rt)
fig

# Three macroscopic attractors: the limit cycle around the origin (blue) and
# the two stable fixed points (red, green). The thin fingers extending outward
# are the stable manifolds of the saddles.

# ## Discretisation: the diffusion generator
#
# Build a `DiffusionGenerator` on the same grid. `DiffusionGenerator` uses the
# Scharfetter–Gummel exponential-fitting finite-volume scheme; with reflecting
# boundaries the discrete generator preserves probability mass.

# The grid eltype is `BigFloat` (default 256-bit). In `Float64` the true
# spectral gap of $Q^\top$ falls below `eps`, so a direct sparse LU
# converges to a quasi-stationary mode rather than the global invariant.
# `SparspakFactorization` is the sparse direct LU that accepts arbitrary-
# precision matrices; swap `BigFloat` for `Double64` (~6× faster) if
# $D \gtrsim 10^{-3}$ is enough for your problem.

grid = CartesianGrid{BigFloat}((-5.5, 5.5, 301), (-3.0, 3.0, 301))
gen = DiffusionGenerator(sys, grid)
@show size(gen.Q)
@show Float64(maximum(abs, vec(sum(gen.Q; dims = 2))))  # row sums ≈ 0

# Cell centers as `Float64` for plotting; the underlying solve stays in
# `BigFloat`.
reshape_grid(v) = reshape(v, grid.nbox)
xs_c = Float64.(grid.centers[1])
ps_c = Float64.(grid.centers[2])

# ## Quasipotential
#
# $\Phi(x, p) = -\lim_{D\to 0} D\log\rho_{ss}^D$. At finite $D$ we plot
# $\Phi \simeq -D\log\rho_{ss}$ directly. The fixed-point wells dominate;
# the limit cycle is metastable but sits as a *local* min of the global
# $\Phi$ at height $\approx \Delta\Phi_\text{FP→saddle}$ — invisible on a
# full-range colorbar but easy to see when zoomed (see below).

ρ = stationary_distribution(gen, LS.SparspakFactorization())

Φ = .-(D .* log.(max.(ρ, eps(eltype(ρ)))))
Φ .-= minimum(Φ)
@show Float64.(extrema(Φ))

Φg = Float64.(reshape_grid(Φ))
ρg = Float64.(reshape_grid(ρ))
qhi = quantile(filter(isfinite, vec(Φg)), 0.999)
ρhi = quantile(vec(ρg), 0.999)

fig = Figure(; size = (900, 900))

ax_ρ = Axis(
    fig[1, 1]; xlabel = "x", ylabel = "p", aspect = DataAspect(),
    title = "Stationary density ρ_ss  (D = $D)"
)
hm_ρ = heatmap!(ax_ρ, xs_c, ps_c, ρg; colormap = :magma, colorrange = (0, ρhi))
scatter!(ax_ρ, [-x_fp, x_fp], [0.0, 0.0]; color = :white, markersize = 10)
scatter!(
    ax_ρ, [-x_saddle, x_saddle], [0.0, 0.0]; color = :orange,
    marker = :diamond, markersize = 8
)
Colorbar(fig[1, 2], hm_ρ; label = "ρ_ss")

ax_Φ = Axis(
    fig[2, 1]; xlabel = "x", ylabel = "p", aspect = DataAspect(),
    title = "Quasipotential Φ = -D log ρ_ss  (D = $D)"
)
hm_Φ = heatmap!(ax_Φ, xs_c, ps_c, Φg; colormap = :viridis, colorrange = (0, qhi))
contour!(ax_Φ, xs_c, ps_c, Φg; levels = 14, color = (:white, 0.4), linewidth = 0.7)
# streamplot!(
#     ax_Φ, u -> Point2f(u[2], -u[2] * 𝒟(u[1]) - V′(u[1])),
#     -5.5 .. 5.5, -3.0 .. 3.0; gridsize = (28, 14), arrow_size = 6,
#     stepsize = 0.005, linewidth = 0.5,
#     colormap = [(:gray85, 0.6), (:gray85, 0.6)]
# )
scatter!(ax_Φ, [-x_fp, x_fp], [0.0, 0.0]; color = :red, markersize = 12)
scatter!(
    ax_Φ, [-x_saddle, x_saddle], [0.0, 0.0]; color = :orange,
    marker = :diamond, markersize = 10
)
Colorbar(fig[2, 2], hm_Φ; label = "Φ")
fig

# **Left:** the raw stationary density. On a linear scale ρ_ss is a pair of
# delta-like spikes localised on the two fixed-point wells — essentially all
# the equilibrium mass sits within `r ~ √D` of $(\pm\sqrt{20}, 0)$, and the
# rest of phase space (limit cycle included) is indistinguishable from
# zero. **Right:** the same data viewed as Φ = -D log ρ_ss. The two fixed
# points sit at the bottom of deep wells; the saddles at $(\pm\sqrt{5}, 0)$
# are the lowest points on the ridge between the basins, which the optimal
# (instanton) escape path crosses. The limit-cycle well is very shallow on
# this scale — see the committor plot below for a finer visualisation of
# its role.

# ## Instanton limit ($D \to 0$)
#
# Everything above is a *finite-D* calculation: the PDE solves for $\rho_{ss}$
# at the actual $D$ we picked, and $\Phi = -D \log \rho_{ss}$ inherits
# $O(D)$ corrections to its strict $D \to 0$ value. The Freidlin–Wentzell
# theory delivers that limit directly: the *anchored* quasipotential
# $\Phi(x;\, A) = \inf_\gamma S[\gamma]$ over paths $\gamma$ from attractor
# $A$ to $x$ is the minimum action of the instanton, computed pointwise by
# the geometric minimum-action method (gMAM). gMAM is $D$-independent;
# repeating it for many endpoints sweeps out $\Phi(x; A)$ over phase space.
#
# We anchor at the left fixed point FP⁻ = $(-\sqrt{20}, 0)$ and sweep a
# coarse sub-grid over the same $(x, p)$ box used for the PDE. Each
# endpoint requires its own minimisation, so the runtime scales linearly
# with the number of sub-grid points.
#
# Two caveats. **(1)** The Keldysh noise covariance $\Sigma =
# \mathrm{Diag}(0, 2D)$ is rank-deficient (no x-channel noise); we
# regularise with a tiny $\epsilon$ on the x-channel so
# $\Sigma_\epsilon = \mathrm{Diag}(\epsilon^2, 2D)$ becomes invertible.
# **(2)** gMAM returns the *geometric* action $S_\text{geo} = 2\Phi$, so we
# divide by 2 throughout.

ε_reg = 0.01
diff_reg(u, p, t) = SA[ε_reg     0.0;
                       0.0   sqrt(2 * D)]
sys_reg = CoupledSDEs(drift_kel, [-0.5, 0.0]; g = diff_reg, noise_strength = 1.0)
x_fp_neg = -x_fp

# Sweep a coarse sub-grid spanning the same phase-space box, computing the
# gMAM action from FP⁻ to each sub-grid point. The result is the
# *anchored* quasipotential $\Phi(x; \text{FP}^-)$, which equals the
# global $\Phi$ on the FP⁻ basin and exceeds it on the FP⁺ side
# (where the path must cross the saddle). The sub-grid is coarse because
# each endpoint requires its own gMAM solve; the cost scales linearly with
# the number of sub-grid points.

nx_sub, np_sub = 101, 51
xs_sub = range(-5.0, 5.0; length = nx_sub)
ps_sub = range(-2.5, 2.5; length = np_sub)

diff_r_sweep(u, p, t) = SA[ε_reg   0.0;
                            0.0    sqrt(2 * D)]
sys_r_sweep = CoupledSDEs(drift_kel, [-0.5, 0.0]; g = diff_r_sweep, noise_strength = 1.0)

function _crude_init(x_end, p_end; npoints = 80)
    Δx = x_end - x_fp_neg
    sg = range(0, 1; length = npoints)
    base = Δx == 0 ? 1.0 : Δx
    return reduce(
        hcat, [SA[x_fp_neg + Δx * s, s * p_end + base * sin(π * s)] for s in sg],
    )
end

# Linear-shift warp: parent's converged path → init for a nearby endpoint.
function _warp_path(path::Matrix{Float64}, new_end)
    n = size(path, 2)
    shift = collect(new_end) .- path[:, end]
    out = copy(path)
    @inbounds for k in 1:n
        s = (k - 1) / (n - 1)
        out[:, k] .+= s .* shift
    end
    return out
end

function _gmam_warm(init_path; reltol = 1.0e-6)
    gmr = minimize_geometric_action(
        sys_r_sweep, init_path; reltol = reltol, maxiters = 10_000, show_progress = false,
    )
    return gmr.action / 2, Matrix(Matrix(gmr.path)')
end

# Sweep endpoints in order of increasing distance from FP⁻; each gMAM solve
# warm-starts from the previously-converged nearest neighbour. This propagates
# the path topology smoothly outward instead of relying on a crude init for
# every endpoint.
all_ij = vec([(i, j) for i in 1:nx_sub, j in 1:np_sub])
dists = [hypot(xs_sub[i] - x_fp_neg, ps_sub[j]) for (i, j) in all_ij]
ij_sorted = all_ij[sortperm(dists)]

function _nearest_visited(i, j, paths, Φ; max_action = 6.0)
    # Skip parents whose own gMAM converged to a bad action — those paths
    # carry a topology error that propagates to descendants. Falling back
    # to a crude init is preferable.
    best_idx = nothing; best_d = Inf
    @inbounds for jj in 1:size(paths, 2), ii in 1:size(paths, 1)
        paths[ii, jj] === nothing && continue
        (isnan(Φ[ii, jj]) || Φ[ii, jj] > max_action) && continue
        d = (ii - i)^2 + (jj - j)^2
        if d < best_d; best_d = d; best_idx = (ii, jj); end
    end
    return best_idx
end

println("Sweeping gMAM over $(nx_sub)×$(np_sub) = $(nx_sub * np_sub) sub-grid points "
      * "(sequential, warm-started)...")
Φ_gmam = fill(NaN, nx_sub, np_sub)
paths_gmam = Array{Union{Nothing, Matrix{Float64}}}(undef, nx_sub, np_sub)
fill!(paths_gmam, nothing)

@showprogress for (i, j) in ij_sorted
    x_end, p_end = xs_sub[i], ps_sub[j]
    if hypot(x_end - x_fp_neg, p_end) < 0.05    # endpoint = FP⁻
        Φ_gmam[i, j] = 0.0
        paths_gmam[i, j] = _crude_init(x_end, p_end)
        continue
    end
    parent = _nearest_visited(i, j, paths_gmam, Φ_gmam)
    init_p = parent === nothing ? _crude_init(x_end, p_end) :
             _warp_path(paths_gmam[parent[1], parent[2]], [x_end, p_end])
    Φ_gmam[i, j], paths_gmam[i, j] = _gmam_warm(init_p)
end

# True FW max in this box ≈ h₁ + h₂ ≈ 5.5 (FP⁻ → FP⁺ via both saddles).
# Anything well above ~8 is the optimiser failing to find a competitive
# path. Mask as NaN so the heatmap shows a blank region instead of a yellow
# blob saturating the colormap.
# Φ_gmam[Φ_gmam .> 8.0] .= NaN
# println("  range (NaN above 8): ", extrema(filter(!isnan, Φ_gmam)))

# Side-by-side: PDE Φ (full grid) vs gMAM Φ anchored at FP⁻ (sub-grid).
fig = Figure(; size = (900, 900))

qhi_l = quantile(filter(isfinite, vec(Φg)), 0.985)
ax1 = Axis(
    fig[1, 1]; xlabel = "x", ylabel = "p", aspect = DataAspect(),
    title = "PDE Φ (D = $D, BigFloat)"
)
hm1 = heatmap!(ax1, xs_c, ps_c, Φg; colormap = :viridis, colorrange = (0, qhi_l))
contour!(ax1, xs_c, ps_c, Φg; levels = 14, color = (:white, 0.3), linewidth = 0.6)
scatter!(ax1, [-x_fp, x_fp], [0.0, 0.0]; color = :red, markersize = 10)
scatter!(
    ax1, [-x_saddle, x_saddle], [0.0, 0.0]; color = :orange,
    marker = :diamond, markersize = 8
)
Colorbar(fig[1, 2], hm1; label = "Φ")

qhi_r = quantile(vec(Φ_gmam), 0.985)
ax2 = Axis(
    fig[2, 1]; xlabel = "x", ylabel = "p", aspect = DataAspect(),
    title = "gMAM Φ(x; FP⁻), instanton limit (anchored)"
)
hm2 = heatmap!(ax2, xs_sub, ps_sub, Φ_gmam; colormap = :viridis, colorrange = (0, qhi_r))
contour!(ax2, xs_sub, ps_sub, Φ_gmam; levels = 10, color = (:white, 0.3), linewidth = 0.6)
scatter!(ax2, [-x_fp, x_fp], [0.0, 0.0]; color = :red, markersize = 10)
scatter!(
    ax2, [-x_saddle, x_saddle], [0.0, 0.0]; color = :orange,
    marker = :diamond, markersize = 8
)
Colorbar(fig[2, 2], hm2; label = "Φ(x; FP⁻)")
fig

# The two heatmaps agree on the **FP⁻ side**: FP⁻ is the deepest well in
# both; the contours in the FP⁻ basin track the PDE contours qualitatively.
# They differ on the **FP⁺ side**: PDE Φ collapses to 0 there (global $\Phi$
# is the min over both attractors), while the anchored gMAM Φ keeps climbing
# past saddle⁻ — that's the FW barrier from FP⁻ to FP⁺. To recover the
# global $\Phi$ from gMAM you would also anchor at FP⁺ and take the
# pointwise minimum.
#
# **Quantitative agreement is approximate**: 5-30% on the FP⁻ side near the
# anchor, worse for endpoints far from FP⁻ or in regions where the path
# topology is hard to seed (e.g. far below the x-axis with x > 0). The
# fundamental obstruction is that the bare gMAM machinery assumes invertible
# noise covariance, while the Keldysh $\Sigma = \mathrm{Diag}(0, 2D)$ is
# rank-deficient. The $\epsilon$ regularisation makes the problem solvable
# but introduces a residual penalty $(\dot x - p)^2 / \epsilon^2$ that the
# optimiser cannot drive to zero. The shape of the landscape is correct;
# the absolute values are an upper bound on the true FW action. A
# proper rank-deficient FW solver (custom Hamiltonian via
# `ExtendedPhaseSpace`, with `pinv(\Sigma)` and a velocity-projection step)
# would close this gap.

# ## Transition Path Theory ($\varepsilon > 0$)
#
# We consider transitions between the two fixed-point attractors:
# $A = \{(x, p) : \|(x, p) - (-\sqrt{20}, 0)\| < r\}$ and similarly for $B$
# centred at $(+\sqrt{20}, 0)$. The limit cycle lies between them and plays the
# role of an obstacle that reactive trajectories must navigate around.

r_attr = 0.4
near_Lfp(u) = norm(u .- SA[-x_fp, 0.0]) < r_attr
near_Rfp(u) = norm(u .- SA[x_fp, 0.0]) < r_attr

result = ReactiveTransition(gen, near_Lfp, near_Rfp; alg = LS.SparspakFactorization())

qplus = forward_committor(result)
qminus = backward_committor(result)

# ### Forward committor
#
# $q^{+}(x, p)$ is the probability of reaching $B$ before $A$. The
# 0.5-isocontour traces the stochastic separatrix.

q2d = Float64.(reshape_grid(qplus))

fig = Figure(; size = (820, 460))
ax = Axis(
    fig[1, 1]; xlabel = "x", ylabel = "p", aspect = DataAspect(),
    title = "Forward committor q⁺"
)
hm = heatmap!(ax, xs_c, ps_c, q2d; colormap = :coolwarm)
contour!(ax, xs_c, ps_c, q2d; levels = [0.5], color = :black, linewidth = 2)
scatter!(ax, [-x_fp, x_fp], [0.0, 0.0]; color = :black, markersize = 12)
Colorbar(fig[1, 2], hm; label = "q⁺")
fig

# The spiral structure inside the limit-cycle region is striking: the
# committor varies smoothly there because the deterministic dynamics wind
# around the limit cycle, sampling both saddles roughly equally on the way
# out. The separatrix passes through the unstable origin.

# ### Reactive current
#
# The reactive current $J_R$ points from $A$ to $B$ everywhere and is
# divergence-free except on $A\cup B$. The node-centered form returned by
# `reactive_current` is convenient for an arrow plot.

J_nodes, _ = reactive_current(result)
Jx, Jy = J_nodes

step = 8
Is = 1:step:grid.nbox[1]
Js = 1:step:grid.nbox[2]
xq = [xs_c[i] for i in Is, _ in Js]
yq = [ps_c[j] for _ in Is, j in Js]
ux = [Float64(Jx[i, j]) for i in Is, j in Js]
uy = [Float64(Jy[i, j]) for i in Is, j in Js]
mag = sqrt.(ux .^ 2 .+ uy .^ 2)
mmax = maximum(mag)

fig = Figure(; size = (820, 460))
ax = Axis(
    fig[1, 1]; xlabel = "x", ylabel = "p", aspect = DataAspect(),
    title = "Reactive current J_R"
)
heatmap!(ax, xs_c, ps_c, q2d; colormap = :coolwarm, alpha = 0.7)
arrows2d!(
    ax, vec(xq), vec(yq),
    vec(ux ./ (mmax + eps())), vec(uy ./ (mmax + eps()));
    color = :black, lengthscale = 0.18
)
fig

# ### Mean first-passage time
#
# Time to reach either fixed-point attractor as a function of starting
# position. The limit-cycle interior is the region of longest wait.

target = u -> near_Lfp(u) || near_Rfp(u)
τ = mean_first_passage_time(gen, target; alg = LS.SparspakFactorization())
τ2d = Float64.(reshape_grid(τ))

fig = Figure(; size = (820, 460))
ax = Axis(
    fig[1, 1]; xlabel = "x", ylabel = "p", aspect = DataAspect(),
    title = "Mean first-passage time to either FP"
)
hm = heatmap!(ax, xs_c, ps_c, τ2d; colormap = :plasma)
contour!(ax, xs_c, ps_c, τ2d; levels = 8, color = (:white, 0.3), linewidth = 0.7)
scatter!(
    ax, [p[1] for p in lc_pts], [p[2] for p in lc_pts];
    color = :white, markersize = 1.5
)
Colorbar(fig[1, 2], hm; label = "τ")
fig

# Inside the limit-cycle basin τ saturates at a plateau set by the escape
# rate over the saddle barrier; on the basins of the fixed points the system
# falls in within a few oscillation periods.

# ## Reusing the generator
#
# Because `DiffusionGenerator` is decoupled from the metastable sets, the same
# sparse matrix can be queried for many different $A, B$ pairs at the cost of
# one linear solve each. For example, sweeping the FP-attractor cutoff radius:

for r in (0.3, 0.4, 0.5, 0.6)
    local Ar = u -> norm(u .- SA[-x_fp, 0.0]) < r
    local Br = u -> norm(u .- SA[x_fp, 0.0]) < r
    local res_r = ReactiveTransition(gen, Ar, Br; alg = LS.SparspakFactorization())
    println("r = $r → reactive_rate ≈ $(Float64(reactive_rate(res_r)))")
end

# At this small noise the absolute reactive rate underflows (Kramers
# escape factor $e^{-\Delta\Phi/D}$ with $\Delta\Phi/D \approx 42$), but the
# *shape* of the committor and the reactive current are well-defined and
# numerically stable — they characterise the geometry of the transition, not
# its absolute rate.
