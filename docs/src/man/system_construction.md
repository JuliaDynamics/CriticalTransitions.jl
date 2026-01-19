# Defining a forced dynamical system

This section explains how to specify your dynamical system and forcing of interest.

To specify a system, `CriticalTransitions` provides three core system types:

- [CoupledODEs](@ref) - used to define a deterministic system of ordinary differential equations 
  of the form ``\frac{\text{d}\mathbf{u}}{\text{d}t} = \mathbf{f}(\mathbf{u}, p, t)``.
- [CoupledSDEs](@ref) - used to define a system of stochastic differential equations
  of the form ``\text{d}\mathbf{u} = \mathbf{f}(\mathbf{u}, p, t) \text{d}t + \mathbf{g}(\mathbf{u}, p, t) \text{d}\mathcal{N}_t``.
- [RateSystem](@ref) - used to define a non-autonomous system with parametric forcing
  of the form ``\frac{\text{d}\mathbf{u}}{\text{d}t} = \mathbf{f}(\mathbf{u}(t), p(t))``.

The `CoupledODEs` and `CoupledSDEs` system types are inherited from [`DynamicalSystemsBase.jl`](https://juliadynamics.github.io/DynamicalSystemsDocs.jl/dynamicalsystemsbase/stable/).
The `RateSystem` type is added in CriticalTransitions.jl to enable easy construction of non-autonomous dynamical systems in which a parameter changes over time.

## Deterministic: `CoupledODEs`

```@docs
CoupledODEs
```

## Stochastic: `CoupledSDEs`

```@docs
CoupledSDEs
```

!!! tip
    Check out some examples of how to construct different types of stochastic dynamics [here](@ref defining-stochastic-dynamics).

!!! info
    Note that nonlinear mixings of the Noise Process $\mathcal{W}$ fall into the class of random ordinary differential equations (RODEs) which have a separate set of solvers. See [this example](https://docs.sciml.ai/DiffEqDocs/stable/tutorials/rode_example/) of DifferentialEquations.jl.

### `CoupledSDEs` API

```@docs
solver
drift
div_drift
StochasticSystemsBase.covariance_matrix
StochasticSystemsBase.diffusion_matrix
noise_process
```

## Non-autonomous: `RateSystem`

```@docs
RateSystem
ForcingProfile
frozen_system
```

![Schematic explaining RateSystem construction](../assets/ratesystem_scheme.png)

### `RateSystem` API

```@docs
parameters
set_forcing_start!
set_forcing_duration!
set_forcing_scale!
```

## Converting between systems

The deterministic part of a [`CoupledSDEs`](@ref) can be extracted as a 
[`CoupledODEs`](@ref), making it compatible with functionality of `DynamicalSystems.jl`.
In turn, a `CoupledODEs` can easily be extended into a `CoupledSDEs`.

```@docs
CoupledODEs(ds::CoupledSDEs; kwargs...)
CoupledSDEs(ds::CoupledODEs, p; kwargs...)
```

For example, the
Lyapunov spectrum of a `CoupledSDEs` in the absence of noise, here exemplified by the
FitzHugh-Nagumo model, can be computed by typing:

```@example
using CriticalTransitions
using ChaosTools: lyapunovspectrum

function fitzhugh_nagumo(u, p, t)
    x, y = u
    ϵ, β = p
    dx = (-x^3 + x - y)/ϵ
    dy = -β*y + x
    return SVector{2}([dx, dy])
end

p = [1.,3.]

# Define the CoupledSDEs
sys = CoupledSDEs(fitzhugh_nagumo, ones(2), p; noise_strength=0.1)

# Compute Lyapunov spectrum of deterministic flow
ls = lyapunovspectrum(CoupledODEs(sys), 10000)
```