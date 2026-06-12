using CriticalTransitions
using StaticArrays

using Logging

# Maier-Stein β = 10: drift is a 2D vector field with three fixed points,
# attractors at (±1, 0) and a saddle at (0, 0). The instanton path from
# attractor → attractor crosses the saddle, so it cannot be solved as a
# single shooting BVP (the arclength reparameterization develops a singularity
# wherever H_p = 0, which on the H = 0 shell forces (x, p) = (x*, 0) at any
# drift fixed point). The user must split into two `attractor → saddle` legs.
function maier_stein(u, p, t)
    x, y = u
    return SA[x - x^3 - 10 * x * y^2, -(1 + x^2) * y]
end

ds = CoupledSDEs(maier_stein, zeros(2); noise_strength = 0.25)
H = FreidlinWentzellHamiltonian(ds)

# Reference: gMAM on the full attractor → attractor path.
Nt = 30
xx = collect(range(-1.0, 1.0; length = Nt))
yy = 0.3 .* (1 .- xx .^ 2)
x_initial_full = Matrix([xx yy]')

res_gmam = minimize_geometric_action(
    H, x_initial_full, GeometricGradient(; stepsize = 1.0);
    maxiters = 500, show_progress = false,
)
println(
    "gMAM action (full path, deterministic return makes second half free):  ",
    res_gmam.action
)

# Shooting: left half (-1, 0) → (0, 0).
xs_left = collect(range(-1.0, 0.0; length = Nt))
ys_left = 0.3 .* sin.(π .* (xs_left .+ 1))
x_init_left = Matrix([xs_left ys_left]')

res_left = with_logger(NullLogger()) do
    minimize_geometric_action(
        H, x_init_left,
        MultipleShooting(; nshoots = 8, maxiters = 200, abstol = 1.0e-6);
        show_progress = false,
    )
end
println("Shooting action (-1, 0) → (0, 0):  ", res_left.action)

# Right half by x → -x symmetry: (+1, 0) → (0, 0) has equal action.
xs_right = collect(range(1.0, 0.0; length = Nt))
ys_right = 0.3 .* sin.(π .* (1 .- xs_right))
x_init_right = Matrix([xs_right ys_right]')

res_right = with_logger(NullLogger()) do
    minimize_geometric_action(
        H, x_init_right,
        MultipleShooting(; nshoots = 8, maxiters = 200, abstol = 1.0e-6);
        show_progress = false,
    )
end
println("Shooting action (+1, 0) → (0, 0):  ", res_right.action)

# In Freidlin-Wentzell theory the saddle → attractor leg follows the
# deterministic drift and contributes zero action. The full transition
# action is therefore the maximum of the two halves (in the symmetric case
# they are equal). gMAM's full-path action above should match either half.
println(
    "Symmetry check: |left - right| / max = ",
    abs(res_left.action - res_right.action) / max(res_left.action, res_right.action)
)
