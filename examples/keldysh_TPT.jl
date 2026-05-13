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

using CriticalTransitions, StaticArrays
import LinearSolve as LS
using Attractors
using LinearAlgebra: norm
using Statistics: quantile
using CairoMakie

const ω₀ = 1.0
const γ = 0.03
const μ = 0.1
const α₃ = -0.25
const α₅ = 0.01
const D = 0.1  # noise strength D in the action; momentum noise = √(2D)

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

grid = CartesianGrid((-5.5, 5.5, 251), (-3.0, 3.0, 251))
gen = DiffusionGenerator(sys, grid)
@show size(gen.Q)
@show maximum(abs, vec(sum(gen.Q; dims = 2)))  # row sums ≈ 0

reshape_grid(v) = reshape(v, grid.nbox)
xs_c = collect(grid.centers[1])
ps_c = collect(grid.centers[2])

# ## Quasipotential ($\varepsilon \to 0$)
#
# The Freidlin–Wentzell quasipotential is
# ```math
# \Phi(x, p) = -\lim_{D\to 0} D\log\rho_{ss}^D(x, p) .
# ```
# At the small but finite $D=10^{-3}$ used here, $\Phi \simeq -D\log\rho_{ss}$
# is dominated by the two deep fixed-point wells: their basins carry
# essentially all of the equilibrium mass and the limit cycle is metastable
# rather than a global minimum of $\Phi$.

ρ = stationary_distribution(gen, LS.UMFPACKFactorization())
Φ = .-(D .* log.(max.(ρ, eps(Float64))))
Φ .-= minimum(Φ)
@show extrema(Φ)

Φg = reshape_grid(Φ)
qhi = quantile(filter(isfinite, vec(Φg)), 0.985)

fig = Figure(; size = (820, 460))
ax = Axis(
    fig[1, 1]; xlabel = "x", ylabel = "p", aspect = DataAspect(),
    title = "Quasipotential Φ = -D log ρ_ss  (D = $D)"
)
hm = heatmap!(ax, xs_c, ps_c, Φg; colormap = :viridis, colorrange = (0, qhi))
contour!(ax, xs_c, ps_c, Φg; levels = 14, color = (:white, 0.4), linewidth = 0.7)
streamplot!(
    ax, u -> Point2f(u[2], -u[2] * 𝒟(u[1]) - V′(u[1])),
    -5.5 .. 5.5, -3.0 .. 3.0; gridsize = (28, 14), arrow_size = 6,
    stepsize = 0.005, linewidth = 0.5,
    colormap = [(:gray85, 0.6), (:gray85, 0.6)]
)
scatter!(ax, [-x_fp, x_fp], [0.0, 0.0]; color = :red, markersize = 12)
scatter!(
    ax, [-x_saddle, x_saddle], [0.0, 0.0]; color = :orange,
    marker = :diamond, markersize = 10
)
Colorbar(fig[1, 2], hm; label = "Φ")
fig

# The two fixed points sit at the bottom of deep wells; the saddles at
# $(\pm\sqrt{5}, 0)$ are the lowest points on the ridge between the basins,
# which the optimal (instanton) escape path crosses. The limit-cycle well is
# very shallow on this scale — see the committor plot below for a finer
# visualisation of its role.

# ## Transition Path Theory ($\varepsilon > 0$)
#
# We consider transitions between the two fixed-point attractors:
# $A = \{(x, p) : \|(x, p) - (-\sqrt{20}, 0)\| < r\}$ and similarly for $B$
# centred at $(+\sqrt{20}, 0)$. The limit cycle lies between them and plays the
# role of an obstacle that reactive trajectories must navigate around.

r_attr = 0.4
near_Lfp(u) = norm(u .- SA[-x_fp, 0.0]) < r_attr
near_Rfp(u) = norm(u .- SA[x_fp, 0.0]) < r_attr

result = ReactiveTransition(gen, near_Lfp, near_Rfp)

qplus = forward_committor(result)
qminus = backward_committor(result)

# ### Forward committor
#
# $q^{+}(x, p)$ is the probability of reaching $B$ before $A$. The
# 0.5-isocontour traces the stochastic separatrix.

q2d = reshape_grid(qplus)

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
ux = [Jx[i, j] for i in Is, j in Js]
uy = [Jy[i, j] for i in Is, j in Js]
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
τ = mean_first_passage_time(gen, target)
τ2d = reshape_grid(τ)

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
    local res_r = ReactiveTransition(gen, Ar, Br)
    println("r = $r → reactive_rate ≈ $(reactive_rate(res_r))")
end

# At this small noise the absolute reactive rate underflows (Kramers
# escape factor $e^{-\Delta\Phi/D}$ with $\Delta\Phi/D \approx 42$), but the
# *shape* of the committor and the reactive current are well-defined and
# numerically stable — they characterise the geometry of the transition, not
# its absolute rate.
