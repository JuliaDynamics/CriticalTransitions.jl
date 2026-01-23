
# # String method on the Muller-Brown landscape using a gradient function

# This example demonstrates how to use `CriticalTransitions.string_method` when you *do not*
# have a `DynamicalSystem`, but instead have access to a potential landscape $V(x)$ and its
# gradient $\nabla V(x)$. Together with the implementation of the string method, this examples is inspired by [`Gideon Simpson`](https://www.math.drexel.edu/~simpson/) package [`StringMethod.jl`](https://github.com/gideonsimpson/StringMethod.jl).
#
# The string method implementation in CriticalTransitions expects a drift (vector field)
# `b(x)`. For a gradient system, the natural choice is the gradient-descent flow
#
# ```math
# \dot{x} = -\nabla V(x),
# ```
#
# so we pass `b(x) = -gradV(x)`.

using CriticalTransitions

using CairoMakie
using ForwardDiff
using LinearAlgebra: norm
using NLsolve
using Optim
using Statistics: mean

# ## Muller-Brown potential

# We use the standard 4-Gaussian Muller-Brown potential
#
# ```math
# V(x,y) = \sum_{i=1}^4 A_i \exp\left(a_i(x-x_i)^2 + b_i(x-x_i)(y-y_i) + c_i(y-y_i)^2\right).
# ```

function Muller(x)
    aa = (-1.0, -1.0, -6.5, 0.7)
    bb = (0.0, 0.0, 11.0, 0.6)
    cc = (-10.0, -10.0, -6.5, 0.7)
    AA = (-200.0, -100.0, -170.0, 15.0)
    XX = (1.0, 0.0, -0.5, -1.0)
    YY = (0.0, 0.5, 1.5, 1.0)

    x1, x2 = x
    V = 0.0
    @inbounds for i in 1:4
        dx = x1 - XX[i]
        dy = x2 - YY[i]
        V += AA[i] * exp(aa[i] * dx^2 + bb[i] * dx * dy + cc[i] * dy^2)
    end
    return V
end

V = x -> Muller(x)

# The gradient can be computed with ForwardDiff.
gradV = x -> ForwardDiff.gradient(V, x)

# ## Locate minima and saddles (for plotting)

# This part is not strictly necessary for computing the string, but it helps for
# visualization and for choosing endpoints.

min1 = optimize(V, [-0.75, 1.5])
x1 = min1.minimizer
min2 = optimize(V, [0.6, 0.0])
x2 = min2.minimizer
min3 = optimize(V, [0.0, 0.5])
x3 = min3.minimizer

saddle1 = nlsolve(gradV, [-1.0, 0.5])
s1 = saddle1.zero
saddle2 = nlsolve(gradV, [0.25, 0.3])
s2 = saddle2.zero

x1, x2, x3, s1, s2

# ## Run the string method using a gradient function

# The string method expects a drift field, so we wrap the gradient into
# gradient-descent drift.
b(x) = -gradV(x)

# Choose endpoints between two minima (left well to right well).
xa = x1
xb = x2

# Construct an initial string as a 2×Nt matrix, with a small transverse perturbation
# to make the convergence visually obvious.
function linear_string(a, b, Nt)
    xs = range(a[1], b[1]; length=Nt)
    ys = range(a[2], b[2]; length=Nt)
    return vcat(xs', ys')
end

Nt = 50
x_initial = linear_string(xa, xb, Nt)
x_initial[2, :] .+= 0.15 .* sin.(range(0, π; length=Nt))

# The parameter `stepsize` is the (fixed) time step used internally for the evolution step.
# The Muller-Brown potential has steep regions, so small `stepsize` is recommended.
stepsize = 1e-4
maxiters = 2500

string = CriticalTransitions.string_method(
    b, x_initial; stepsize, maxiters, show_progress=false
)

# A simple convergence diagnostic: average step along the string.
string_m = Matrix(string.path)
dmean = mean(norm(string_m[i + 1, :] .- string_m[i, :]) for i in 1:(size(string_m, 1) - 1))
dmean

# ## Visualize on the potential landscape

xx = LinRange(-1.5, 1.5, 250)
yy = LinRange(-0.5, 2.0, 250)
V_vals = [V([x, y]) for y in yy, x in xx]
V_clip = min.(V_vals, 500)

fig = Figure(; size=(700, 450), fontsize=13)
ax = Axis(fig[1, 1]; xlabel="x", ylabel="y", aspect=1.2)
contour!(ax, xx, yy, V_clip; levels=range(-150, 500, 35), colormap=:viridis)

lines!(ax, x_initial[1, :], x_initial[2, :]; color=:black, linewidth=2, linestyle=:dash) # Initial string (dashed) and converged string (solid)
lines!(ax, string_m[:, 1], string_m[:, 2]; color=:black, linewidth=3)

scatter!(ax, [x1[1], x2[1], x3[1]], [x1[2], x2[2], x3[2]]; color=:red, markersize=10) # Minima and saddles
scatter!(ax, [s1[1], s2[1]], [s1[2], s2[2]]; color=:green, marker=:diamond, markersize=10)

limits!(ax, extrema(xx)..., extrema(yy)...)
fig
