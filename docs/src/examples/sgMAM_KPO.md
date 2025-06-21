```@meta
EditURL = "../../../examples/sgMAM_KPO.jl"
```

We demonstrate the simple geometric minimum action method (sgMAM) on the Kerr parametric oscillator (KPO). The method computes the optimal path between two attractors in the phase space that minimizes the action of the system. It is a simplification of the geometric minimum action method (gMAM) by avoiding the computation of the second order derivatives of the extended Hamiltonian of the optimisation problem.

````@example sgMAM_KPO
using CriticalTransitions, CairoMakie
````

The KPO equation is a nonlinear ordinary differential equation that describes the response of the nonlinear parametrically driven resonator at its dominant resonant condition. The equation of motion are of the form:
```math
\dot{\mathbf{x}} = \mathbf{f}(\mathbf{x}) + \sigma\mathbf{ξ}(t)
```
where `f` is an autonemous drift function and and we have brownian noise `ξ` with intensity `σ`.

Here we define the define the drift of each seperable variable `u` and `v`. In addition, we hard-code the Jacobian of the drift function.

````@example sgMAM_KPO
const λ = 3 / 1.21 * 2 / 295
const ω0 = 1.000
const ω = 1.000
const γ = 1 / 295
const η = 0
const α = -1

function fu(u, v)
    return (-4 * γ * ω * u - 2 * λ * v - 4 * (ω0 - ω^2) * v - 3 * α * v * (u^2 + v^2)) /
           (8 * ω)
end
function fv(u, v)
    return (-4 * γ * ω * v - 2 * λ * u + 4 * (ω0 - ω^2) * u + 3 * α * u * (u^2 + v^2)) /
           (8 * ω)
end
stream(u, v) = Point2f(fu(u, v), fv(u, v))
dfvdv(u, v) = (-4 * γ * ω + 6 * α * u * v) / (8 * ω)
dfudu(u, v) = (-4 * γ * ω - 6 * α * u * v) / (8 * ω)
dfvdu(u, v) = (-2 * λ + 4 * (ω0 - ω^2) + 9 * α * u^2 + 3 * α * v^2) / (8 * ω)
dfudv(u, v) = (-2 * λ - 4 * (ω0 - ω^2) - 3 * α * u^2 - 9 * α * v^2) / (8 * ω)
````

The optimisation is performed in a doubled phase space, i.e., every variable of the SDE system is considered as a generelised coordinate $\mathbf{x}$ and gets a corresponding generalised momentum $\mathbf{p}$. The makes it that also systems with dissipative flow can be solved. As such, we extend the phase space by defining the hamiltionian
```math
H = \sum_i \frac{p_i^2}{2} + f_i(\mathbf{x})p_i
```
Hence, to use the sgMAM method, we need to define the derivatives of the Hamiltonian with respect to the phase space coordinates and the generalised momentum:

````@example sgMAM_KPO
function H_x(x, p) # ℜ² → ℜ²
    u, v = eachrow(x)
    pu, pv = eachrow(p)

    H_u = @. pu * dfudu(u, v) + pv * dfvdu(u, v)
    H_v = @. pu * dfudv(u, v) + pv * dfvdv(u, v)
    return Matrix([H_u H_v]')
end
function H_p(x, p) # ℜ² → ℜ²
    u, v = eachrow(x)
    pu, pv = eachrow(p)

    H_pu = @. pu + fu(u, v)
    H_pv = @. pv + fv(u, v)
    return Matrix([H_pu H_pv]')
end

sys = SgmamSystem{false,2}(H_x, H_p)
````

We saved this function in the `SgmamSystem` struct. We want to find the optimal path between two attractors in the phase space. We define the initial trajectory as `wiggle` between the two attractors.

````@example sgMAM_KPO
Nt = 500  # number of discrete time steps
s = collect(range(0; stop=1, length=Nt))

xa = [-0.0208, 0.0991]
xb = -xa
xsaddle = [0.0, 0.0]
````

Initial trajectory

````@example sgMAM_KPO
xx = @. (xb[1] - xa[1]) * s + xa[1] + 4 * s * (1 - s) * xsaddle[1]
yy = @. (xb[2] - xa[2]) * s + xa[2] + 4 * s * (1 - s) * xsaddle[2] + 0.01 * sin(2π * s)
x_initial = Matrix([xx yy]')
````

The optimisation is the performed by the `sgmam` function:

````@example sgMAM_KPO
MLP = sgmam(sys, x_initial; iterations=1_000, ϵ=10e2, show_progress=false)
x_min = MLP.path;
nothing #hide
````

The function returns the optimal path `x_min`, the minmal action `S_min`, the Lagrange multipliers `lambda` associated with the optimal path, the optimal generalised momentum `p`, and the time derivative of the optimal path `xdot`. We can plot the initial trajectory and the optimal path:

````@example sgMAM_KPO
fig, ax, _ = lines(
    x_initial[1, :], x_initial[2, :]; label="init", linewidth=3, color=:black
)
lines!(x_min[:, 1], x_min[:, 2]; label="MLP", linewidth=3, color=:red)
streamplot!(
    ax,
    stream,
    (-0.08, 0.08),
    (-0.15, 0.15);
    gridsize=(20, 20),
    arrow_size=10,
    stepsize=0.001,
    colormap=[:gray, :gray],
)
axislegend(ax)
fig
````

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

