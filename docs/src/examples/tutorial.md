# Tutorial

To give you an idea of how our package works, this tutorial provides some example code with explanations.

## Example: FitzHugh-Nagumo model
Let's consider a simple 2-dimensional dynamical system - the *FitzHugh-Nagumo* model:

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

This works exactly as in [DynamicalSystems.jl](https://juliadynamics.github.io/DynamicalSystemsBase.jl/dev/) by defining a function `f(u,p,t)` which takes as input a vector `u` of state variables (``u``,``v``), a vector `p` of parameters, and time `t`. The function must return an array of flow increments ($\text{d}u$, $\text{d}v$). For performance reasons, it is advisable to return a StaticArray `SA[du, dv]` rather than just a Vector `[du, dv]`.

```@example MAIN
using CriticalTransitions
import Random # hide
Random.seed!(1) # hide

function fitzhugh_nagumo(u,p,t)
    u, v = u
    ϵ, β, α, γ, κ, I = p

    du = (-α*u^3 + γ*u - κ*v + I)/ϵ
    dv = -β*v + u

    SA[du, dv]
end
```

!!! tip "In-place vs. out-of-place"
    The function `fitzhugh_nagumo(u,p,t)` is defined *out-of-place*. It is also possible to define the system *in-place* as `fitzhugh_nagumo!(du,u,p,t)`. For more info, see [here](https://diffeq.sciml.ai/stable/types/ode_types/).

### CoupledSDE

Next, we construct a stochastic system with the `fitzhugh_nagumo` equation as the deterministic part. Suppose we would like to force both state variables ``u`` and ``v`` with additive, uncorrelated Gaussian noise of intensity ``\sigma``. This is the default case. We simply write

```@example MAIN
p = [1., 3., 1., 1., 1., 0.] # Parameters (ϵ, β, α, γ, κ, I)
σ = 0.2 # noise strength

# CoupledSDE
sys = CoupledSDEs(fitzhugh_nagumo, zeros(2), p; noise_strength=σ)
```
Here the first field `fitzhugh_nagumo` specifies the deterministic dynamics `f` (see [`CoupledSDEs`](@ref)). We have chosen `zeros(2)` as the initial state of the system, which is the second field. The length of this vector must match the system's dimensionality. In the (optional) third field, we specify the parameter vector `p`, which includes the parameters of `f` followed by the parameters of `g` (in this case, there are no parameters for `g`). Lastly, `noise_strength` sets the noise strength. Since we have not specified a noise process, the default case of an uncorrelated Wiener process is used.

!!! note "Multiplicative and/or correlated noise"
    Of course, it is also possible to define more complicated noise processes than simple additive white noise. This is done by specifying a custom *noise function* and *covariance matrix* in the `CoupledSDEs` definition. For more info, see [`CoupledSDEs`](@ref).

That's it! Now we can apply the toolbox of `CriticalTransitions` to our stochastic FitzHugh-Nagumo system `sys`.

### Find stable equilibria
For the parameters chosen above, the FitzHugh-Nagumo system is bistable. Let's compute the fixed points using the [`ChaosTools.fixedpoints`](@ref) function. This function is borrowed from ChaosTools.jl and is loaded as an extension when we write `using ChaosTools`.

```@example MAIN
using ChaosTools
# Calculate fixed points and store the stable ones
eqs, eigs, stab = fixedpoints(sys, [-2,-2], [2,2])
fp1, fp2 = eqs[stab]
```

### Stochastic simulation
Using the [`trajectory`](@ref) function, we now run a simulation of our system for `200` time units starting out from the fixed point `fp1`:

```@example MAIN
sim = trajectory(sys, 200, fp1)
```

In the keyword arguments, we have specified at which interval the solution is saved. Further keyword arguments can be used to change the solver (the default is `SOSRA()` for stochastic integration) and other settings.

The simulated trajectory is stored in `sim` in the usual output format of the [`solve`](https://docs.sciml.ai/DiffEqDocs/stable/basics/common_solver_opts/) method of DifferentialEquations.jl, including the solution `sim.u` and the vector of time points `sim.t`. The solution can also be accessed as a matrix `sim[i, t]`, where `i` is the `i`-th component of `u` and `t` the time index.

Let's plot the result. Did the trajectory transition to the other attractor?

```@example MAIN
using Plots
plt = plot(sim[1][:,1], sim[1][:,2]; xlabel="u", ylabel="v", legend=false)
scatter!([fp1[1], fp2[1]], [fp1[2], fp2[2]], color=:red, markersize=4)
xlims!(-1.2, 1.2)
ylims!(-0.6, 0.6)
plt
```

Hopefully, this helped you to get started. For more info, check out the Manual section of these docs.