# Tutorial

To give you an idea of how our package works, this tutorial provides some example code with explanations.

## Example: FitzHugh-Nagumo model
Consider the FitzHugh-Nagumo model,

```math
\begin{aligned}
\frac{du}{dt} &= \frac{1}{\epsilon} \left( -\alpha u^3 + \gamma u - \kappa v + I \right) \\
\frac{dv}{dt} &= -\beta v + u \ ,
\end{aligned}
```

where ``\epsilon`` is the parameter of time scale separation between the state variables ``u`` and ``v``. The parameters ``\alpha >0``, ``\beta >1``, ``\gamma>0``, and ``\kappa>0`` are real constants, and ``I`` denotes a driving term.

Let's investigate this system under stochastic forcing.

### System definition
First, we need to translate the system equations above into Julia code.

This works by defining a function `f(u,p,t)` which takes as input a vector `u` of state variables (``u``,``v``), a vector `p` of parameters, and time `t`. The function must return an array of flow increments (``du``, ``dv``). For performance reasons, it is advisable to return a StaticArray `SA[du, dv]` rather than just a Vector `[du, dv]`. This is why we need the `StaticArrays` package.

```@example MAIN
using CriticalTransitions

function fitzhugh_nagumo(u,p,t)
    u, v = u
    ϵ, β, α, γ, κ, I = p[1]

    du = (-α*u^3 + γ*u - κ*v + I)/ϵ
    dv = -β*v + u

    SA[du, dv]
end
```

Note that the system parameters `ϵ, β, α, γ, κ, I = p[1]` are unpacked as the first component of `p`. This is necessary because in CriticalTransitions.jl one can also define a separate set of parameters for the stochastic component of the system, which would then make up the second component `p[2]` ( see [Define a CoupledSDE](@ref)).

!!! tip "In-place vs. out-of-place"
    The function `fitzhugh_nagumo(u,p,t)` is defined *out-of-place*. It is also possible to define the system *in-place* as `fitzhugh_nagumo!(du,u,p,t)`. For more info, see [here](https://diffeq.sciml.ai/stable/types/ode_types/).

### CoupledSDE

Next, we turn the `fitzhugh_nagumo` system into a stochastic dynamical system. Suppose we would like to force both state variables ``u`` and ``v`` with additive, uncorrelated Gaussian noise of intensity ``\sigma``. This is the default case. We simply write

```@example MAIN
p = [1., 3., 1., 1., 1., 0.] # Parameters (ϵ, β, α, γ, κ, I)
σ = 0.18 # noise strength

# CoupledSDE
sys = CoupledSDEs(fitzhugh_nagumo, diag_noise_function(σ), zeros(2), p)
```
Here we have chosen `zeros(2)` as the initial state of the system. The length of this vector must correspond to the system's dimensionality, but for now the state is just a placeholder that aligns our syntax with that of DifferentialEquations.jl and DynamicalSystems.jl.

!!! note "Multiplicative and/or correlated noise"
    Of course, it is also possible to define more complicated noise processes than simple additive white noise. This is done by specifying a custom *noise function* and *covariance matrix* in the `CoupledSDEs` definition. For more info, see [Define a CoupledSDE](@ref).

That's it! Now we can throw the toolbox of `CriticalTransitions` at our stochastic FitzHugh-Nagumo system `sys`.

### Find stable equilibria
For the parameters chosen above, the FitzHugh-Nagumo system is bistable. Let's compute the fixed points using the [`fixedpoints`](https://juliadynamics.github.io/DynamicalSystemsDocs.jl/chaostools/stable/periodicity/#ChaosTools.fixedpoints) function from ChaosTools.jl. As this function is from the `DynamicalSystems` ecosystem, it takes a system of type `CoupledODEs` as input. We can simply convert the CoupledSDEs `sys` via the [`CoupledODEs`](@ref) function:

```@example MAIN
# Calculate fixed points
ds = CoupledODEs(sys)
box = intervals_to_box([-2,-2], [2,2])
eqs, eigs, stab = fixedpoints(ds, box)

# Store the two stable fixed points
fp1, fp2 = eqs[stab]
```

### Stochastic simulation
Using the `simulate` function, we now run a simulation of our system starting out from the fixed point `fp1`:

```@example MAIN
sim = simulate(sys, fp1, dt=0.01, tmax=1e3)
```

In the keyword arguments, we have specified the time step `dt` and total duration `tmax` of the numerical time integration.

The simulated trajectory is stored in `sim` as a matrix with 2 rows corresponding to the state variables ``u``, ``v``, and 10,000 columns corresponding to the time steps.

Let's plot the result. Did the trajectory transition to the other attractor?

```@example MAIN
using CairoMakie
fig = Figure(size = (800, 600), fontsize=18)
ax = Axis(fig[1, 1], xlabel = "u", ylabel = "v")
limits!(-1.2,1.2,-0.6,0.6)
lines!(ax, sim[1,:], sim[2,:])
scatter!(ax, [fp1[1], fp2[1]], [fp1[2], fp2[2]], color=:red, markersize=20)
fig
```

Hopefully, this helped you to get started. For more info, check out the Manual section of these docs.