# Tutorial

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

This works by defining a function `f(u,p,t)` which takes as input a vector `u` of state variables (``u``,``v``), a vector `p` of parameters, and time `t`. The function must return a StaticArray `SA[du, dv]` of flow increments (``du``, ``dv``).

```julia
function fitzhugh_nagumo(u,p,t)
    u, v = u
    ϵ, β, α, γ, κ, I = p[1]

    du = (-α*u^3 + γ*u - κ*v + I)/ϵ
    dv = -β*v + u

    SA[du, dv]
end
```

Returning a StaticArray instead of simply a Vector `[du, dv]` ensures that the above function definition is compatible with DynamicalSystems.jl and improves performance when integrating.

!!! tip "In-place vs. out-of-place"
    The function `fitzhugh_nagumo(u,p,t)` is defined *in-place*. It is also possible to define the system *out-of-place* as `fitzhugh_nagumo!(du,u,p,t)`. For more info, see [here](https://diffeq.sciml.ai/stable/types/ode_types/).

### StochSystem

Next, we turn the `fitzhugh_nagumo` system into a stochastic dynamical system. Suppose we would like to force both state variables ``u`` and ``v`` with additive, uncorrelated Gaussian noise of intensity ``\sigma``. This is the default case. We simply write

```julia
# Parameters (ϵ, β, α, γ, κ, I)
p = [1., 3., 1., 1., 1., 0.]
σ = 0.2

# StochSystem
sys = StochSystem(fitzhugh_nagumo, p, zeros(2), σ)
```
Here we have chosen `zeros(2)` as the initial state of the system. The length of this vector must correspond to the system's dimensionality, but for now the state is just a placeholder that aligns our syntax with that of DifferentialEquations.jl and DynamicalSystems.jl.

!!! note "Multiplicative and/or correlated noise"
    Of course, it is also possible to define more complicated noise processes than simple additive white noise. This is done by specifying a custom *noise function* and *covariance matrix* in the `StochSystem` definition. For more info, see [Defining a StochSystem](@ref).

That's it! Now we can throw the toolbox of `CriticalTransitions` at our stochastic FitzHugh-Nagumo system `sys`.

### Find stable equilibria
For the parameters chosen above, the FitzHugh-Nagumo system is bistable. Let's compute the stable fixed points by calling the [`fixedpoints`](@ref) function.

```julia
# Calculate fixed points
eqs, eigs, stab = fixedpoints(sys, [-2,-2], [2,2])

# Store the two stable fixed points
fp1, fp2 = eqs[stab]
```

### Stochastic simulation
Using the `simulate` function, we now run a simulation of our system starting out from the fixed point `fp1`:

```julia
sim = simulate(sys, fp1, dt=0.01, tmax=1e3)
```

In the keyword arguments, we have specified the time step `dt` and total duration `tmax` of the numerical time integration.

The simulated trajectory is stored in `sim` as a matrix with 2 rows corresponding to the state variables ``u``, ``v``, and 10,000 columns corresponding to the time steps.

Let's plot the result. Did the trajectory tip to the other basin of attraction?

```julia
using PyPlot
scatter([fp1[0], fp2[0]], [fp1[1], fp2[1]], c=["r"])
plot(sim[0,:], sim[1,:])
```

Hopefully, this helped you to get started. For more info, check out the Manual section of these docs.