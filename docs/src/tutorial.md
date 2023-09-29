# Tutorial

To give you an idea of how our package works, this tutorial provides some example code with explanations.

You can just copy and paste the code blocks below, after loading the required packages:

```julia
using CriticalTransitions, StaticArrays
```

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

This works by defining a function `f(u,p,t)` which takes as input a vector `u` of state variables (``u``,``v``), a vector `p` of parameters, and time `t`. The function must return an array of flow increments (``du``, ``dv``). For performance reasons, it is advisable to return a StaticArray `SA[du, dv]` rather than just a Vector `[du, dv]`. This is why we needed the `StaticArrays` package above.

```julia
function fitzhugh_nagumo(u,p,t)
    u, v = u
    ϵ, β, α, γ, κ, I = p[1]

    du = (-α*u^3 + γ*u - κ*v + I)/ϵ
    dv = -β*v + u

    SA[du, dv]
end
```

Note that the system parameters `ϵ, β, α, γ, κ, I = p[1]` are unpacked as the first component of `p`. This is necessary because in CriticalTransitions.jl one can also define a separate set of parameters for the stochastic component of the system, which would then make up the second component `p[2]` ( see [Defining a StochSystem](@ref)).

!!! tip "In-place vs. out-of-place"
    The function `fitzhugh_nagumo(u,p,t)` is defined *out-of-place*. It is also possible to define the system *in-place* as `fitzhugh_nagumo!(du,u,p,t)`. For more info, see [here](https://diffeq.sciml.ai/stable/types/ode_types/).

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
For the parameters chosen above, the FitzHugh-Nagumo system is bistable. Let's compute the fixed points by calling the [`fixedpoints`](@ref) function.

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

Let's plot the result. Did the trajectory transition to the other attractor?

```julia
using PyPlot
scatter([fp1[0], fp2[0]], [fp1[1], fp2[1]], c=["r"])
plot(sim[0,:], sim[1,:])
```

Hopefully, this helped you to get started. For more info, check out the Manual section of these docs.