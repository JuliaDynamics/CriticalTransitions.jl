# # [Transition Path Theory for a 2D double well](@id TPT_example)

# In this example we apply Transition Path Theory (TPT) to a two-dimensional
# overdamped double-well system. The full theoretical background and the layered
# API used below are documented in the
# [Transition Path Theory manual](@ref Transition-Path-Theory). Here we focus on
# putting the pieces together end-to-end: discretise the diffusion on a Cartesian
# grid, solve for the committor, and visualise the reactive density and current.

using CriticalTransitions
using LinearAlgebra: norm
using CairoMakie

# ## System

# We consider an overdamped 2D Itô diffusion
#
# ```math
# \mathrm{d}X_t = -\nabla U(X_t)\,\mathrm{d}t + \sigma\,\mathrm{d}W_t
# ```
#
# with the separable potential
#
# ```math
# U(x, y) = \tfrac{1}{4}(x^2 - 1)^2 + \tfrac{1}{2}\,y^2 ,
# ```
#
# whose drift is
#
# ```math
# b(x, y) = -\nabla U = \bigl(x - x^3,\; -y\bigr) .
# ```
#
# The deterministic dynamics has stable fixed points at $(\pm 1, 0)$ and a saddle
# at the origin. The noise is isotropic with diffusion coefficient $\sigma^2$
# along each axis. Because the drift is the gradient of a potential and the
# noise is isotropic, this is a **reversible** system — the backward committor
# satisfies $q^{-} = 1 - q^{+}$ exactly.

drift(u, p, t) = SA[u[1] - u[1]^3, -u[2]]
sys = CoupledSDEs(drift, [0.0, 0.0]; noise_strength=0.4)

# Visualise the drift field:

xs = range(-1.8, 1.8; length=40)
ys = range(-1.0, 1.0; length=24)
fig = Figure(; size=(640, 360))
ax = Axis(fig[1, 1]; xlabel="x", ylabel="y", aspect=DataAspect())
streamplot!(
    ax,
    p -> Point2f(p[1] - p[1]^3, -p[2]),
    -1.8 .. 1.8,
    -1.0 .. 1.0;
    linewidth=0.6,
    colormap=[:gray40, :gray40],
    arrow_size=8,
    gridsize=(28, 16),
)
scatter!(ax, [-1.0, 1.0], [0.0, 0.0]; color=:black, markersize=12)
text!(ax, [-1.0, 1.0], [0.0, 0.0]; text=["A", "B"], offset=(8, 8))
fig

# ## Discretisation

# We discretise on a uniform Cartesian grid covering the relevant transition
# region. The grid does not have to be square — here we use a finer resolution
# along $x$ since that is the reaction coordinate.

grid = CartesianGrid((-1.8, 1.8, 121), (-1.0, 1.0, 81))

# The metastable sets $A$ and $B$ are passed as predicates evaluated at each
# cell center; alternatively a `BitVector` or a vector of linear cell indices
# is accepted.

A = x -> norm(x .- SA[-1.0, 0.0]) < 0.25
B = x -> norm(x .- SA[1.0, 0.0]) < 0.25

# Building the [`DiffusionGenerator`](@ref) is the expensive step (it sweeps
# the grid once and assembles the sparse rate matrix via the Scharfetter–Gummel
# scheme). Once it exists you can run any number of analyses on the same
# operator.

gen = DiffusionGenerator(sys, grid)

# Sanity-check the discrete generator: rows of a CTMC rate matrix sum to zero
# and the off-diagonals are non-negative.

@show size(gen.Q)
@show maximum(abs, vec(sum(gen.Q; dims=2)))

# ## Bundled TPT solve

# [`ReactiveTransition`](@ref) computes the invariant density, both committors
# and the reactive rate in a single call:

result = ReactiveTransition(gen, A, B)

@show reactive_rate(result)
@show probability_reactive(result)
@show probability_last_A(result)

# `forward_committor`, `backward_committor`, and `stationary_distribution` are
# constant-time accessors on a `ReactiveTransition`:

qplus = forward_committor(result)
qminus = backward_committor(result)
ρ = stationary_distribution(result)

# Reshape from linear indexing into the grid shape for plotting:

reshape_grid(v) = reshape(v, grid.nbox)

# ## Forward committor

# $q^{+}(x, y)$ is the probability that a trajectory started at $(x, y)$
# reaches $B$ before $A$. It interpolates smoothly between $0$ on $A$ and $1$
# on $B$, with the $0.5$-isoline running through the saddle.

x_centers = collect(grid.centers[1])
y_centers = collect(grid.centers[2])

fig = Figure(; size=(640, 360))
ax = Axis(fig[1, 1]; xlabel="x", ylabel="y", aspect=DataAspect(), title="Forward committor q⁺")
hm = heatmap!(ax, x_centers, y_centers, reshape_grid(qplus); colormap=:viridis)
contour!(
    ax, x_centers, y_centers, reshape_grid(qplus); levels=[0.5], color=:white, linewidth=2
)
Colorbar(fig[1, 2], hm)
fig

# Because the system is reversible, $q^{-} = 1 - q^{+}$ — let's verify:

@show maximum(abs.(qminus .- (1 .- qplus)))

# ## Reactive density

# The probability density of being on a reactive piece (one that just left $A$
# and will reach $B$ before returning) is $\rho_R = \rho \cdot q^{+} \cdot q^{-}$.
# It concentrates around the most likely transition channels.

ρR = reactive_density(result)

fig = Figure(; size=(640, 360))
ax = Axis(fig[1, 1]; xlabel="x", ylabel="y", aspect=DataAspect(), title="Reactive density ρ·q⁺·q⁻")
hm = heatmap!(ax, x_centers, y_centers, reshape_grid(ρR); colormap=:magma)
Colorbar(fig[1, 2], hm)
fig

# ## Reactive current

# `reactive_current` returns two `NTuple{D}` outputs: a per-axis array of face
# fluxes (one per axis-aligned face) and a per-axis array of node-centered
# fluxes (averaged onto cell centers). The node-centered version is convenient
# for plotting an arrow field.

J_nodes, _ = reactive_current(result)
Jx, Jy = J_nodes      # each is a 121×81 array

# Subsample for a clean arrow plot:

step = 6
Is = 1:step:grid.nbox[1]
Js = 1:step:grid.nbox[2]
xq = [x_centers[i] for i in Is, _ in Js]
yq = [y_centers[j] for _ in Is, j in Js]
ux = [Jx[i, j] for i in Is, j in Js]
uy = [Jy[i, j] for i in Is, j in Js]
mag = sqrt.(ux .^ 2 .+ uy .^ 2)
mmax = maximum(mag)

fig = Figure(; size=(640, 360))
ax = Axis(fig[1, 1]; xlabel="x", ylabel="y", aspect=DataAspect(), title="Reactive current")
heatmap!(
    ax, x_centers, y_centers, reshape_grid(ρR); colormap=:magma, alpha=0.6,
)
arrows2d!(
    ax,
    vec(xq), vec(yq),
    vec(ux ./ mmax), vec(uy ./ mmax);
    color=:white,
    lengthscale=0.08,
)
fig

# The reactive current points consistently from $A$ toward $B$ and is strongest
# near the saddle, as expected.

# ## Mean first-passage time

# The same [`DiffusionGenerator`](@ref) supports other backward-Kolmogorov
# analyses without rebuilding. For instance, the mean first-passage time to
# $B$ as a function of starting point:

τ = mean_first_passage_time(gen, B)

fig = Figure(; size=(640, 360))
ax = Axis(fig[1, 1]; xlabel="x", ylabel="y", aspect=DataAspect(), title="Mean first-passage time to B")
hm = heatmap!(ax, x_centers, y_centers, reshape_grid(τ); colormap=:plasma)
Colorbar(fig[1, 2], hm)
fig

# ## Reusing the generator

# Because the generator is decoupled from the metastable sets, we can sweep
# over different choices of $A$ and $B$ at the cost of a single sparse solve
# each:

for r in (0.20, 0.25, 0.30, 0.35)
    Ar = x -> norm(x .- SA[-1.0, 0.0]) < r
    Br = x -> norm(x .- SA[1.0, 0.0]) < r
    k = reactive_rate(ReactiveTransition(gen, Ar, Br))
    println("r = $r → k_AB = $(round(k; sigdigits=4))")
end

# Authored by O. Ameye and R. Börner
