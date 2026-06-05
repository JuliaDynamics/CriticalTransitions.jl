# Defining a forced dynamical system

This section explains how to specify your dynamical system and forcing of interest.

To specify a system, `CriticalTransitions` provides three core system types:

- [CoupledODEs](@ref) - used to define a deterministic system of ordinary differential equations of the form ``\frac{\text{d}\mathbf{u}}{\text{d}t} = \mathbf{f}(\mathbf{u}, p, t)``.
- [CoupledSDEs](@ref) - used to define a system of stochastic differential equations of the form ``\text{d}\mathbf{u} = \mathbf{f}(\mathbf{u}, p, t) \text{d}t + \mathbf{g}(\mathbf{u}, p, t) \text{d}\mathcal{N}_t``.
- [RateSystem](@ref) - used to define a non-autonomous system with parametric forcing of the form ``\frac{\text{d}\mathbf{u}}{\text{d}t} = \mathbf{f}(\mathbf{u}(t), p(t))``.

The `CoupledODEs` and `CoupledSDEs` system types are inherited from [`DynamicalSystemsBase.jl`](https://juliadynamics.github.io/DynamicalSystemsDocs.jl/dynamicalsystemsbase/stable/). The `RateSystem` type is added in CriticalTransitions.jl to enable easy construction of non-autonomous dynamical systems in which a parameter changes over time.

## Deterministic: `CoupledODEs`

```@docs
CoupledODEs
```

## Stochastic: `CoupledSDEs`

```@docs; canonical=false
CoupledSDEs
```

!!! tip
    Check out some examples of how to construct different types of stochastic dynamics [here](@ref defining-stochastic-dynamics).

!!! info
    Note that nonlinear mixings of the Noise Process $\mathcal{W}$ fall into the class of random ordinary differential equations (RODEs) which have a separate set of solvers. See [this example](https://docs.sciml.ai/DiffEqDocs/stable/tutorials/rode_example/) of DifferentialEquations.jl.

### `CoupledSDEs` API

```@docs; canonical=false
solver
drift
div_drift
covariance_matrix
diffusion_matrix
noise_process
noise_strength
```

### Noise strength and the trace convention

A [`CoupledSDEs`](@ref) exposes its noise via the diffusion function ``g`` and, for additive invertible noise, the noise rate covariance
```math
\mathbf{\Sigma} \;=\; \texttt{covariance\_matrix(sys)} \;=\; g\,g^\top.
```
For an SDE built as ``\mathrm{d}\mathbf{x} = \mathbf{b}\,\mathrm{d}t + \sigma\,\mathbf{\Sigma}_0\,\mathrm{d}\mathbf{W}_t`` with ``\mathbf{\Sigma}_0\mathbf{\Sigma}_0^\top = \mathbf{Q}``, this gives ``\mathbf{\Sigma} = \sigma^2\,\mathbf{Q}``.

The noise strength and the covariance are not separately determined by the SDE: only their product ``\sigma^2\mathbf{Q}`` enters, and ``(\sigma, \mathbf{Q}) \leftrightarrow (c\sigma, \mathbf{Q}/c^2)`` describes the same physical system. The package picks a unique ``(\sigma, \mathbf{Q})`` pair from ``\mathbf{\Sigma}`` by adopting the **trace convention** ``\mathrm{tr}(\mathbf{Q}) = D``:
```math
\sigma_{\mathrm{eff}}^2 \;=\; \frac{\mathrm{tr}(\mathbf{\Sigma})}{D},
\qquad
\mathbf{Q}_{\mathrm{can}} \;=\; \frac{D}{\mathrm{tr}(\mathbf{\Sigma})}\,\mathbf{\Sigma}.
```
The accessor [`noise_strength`](@ref)`(sys)` returns ``\sigma_{\mathrm{eff}}``. The trace is invariant under orthogonal changes of basis and reduces to the per-direction average noise variance, making this the natural choice among rotation-invariant normalizations.

This convention is what the [large deviation action functionals](@ref "The action") use: they evaluate the Freidlin-Wentzell integrand in the ``\mathbf{Q}_{\mathrm{can}}^{-1}`` metric,
```math
\texttt{fw\_action(sys, path, time)} \;=\; \Phi_{FW}^{\mathbf{Q}_{\mathrm{can}}}[\phi] \;=\;
    \tfrac{1}{2}\!\int (\dot\phi - \mathbf{b})^\top \mathbf{Q}_{\mathrm{can}}^{-1}(\dot\phi - \mathbf{b})\,\mathrm{d}t,
```
so the returned action is independent of the `noise_strength` keyword chosen at construction and invariant under orthogonal changes of basis.

If you built your system as ``\mathrm{d}\mathbf{x} = \mathbf{b}\,\mathrm{d}t + \sigma_{\mathrm{user}}\,\mathbf{\Sigma}_0\,\mathrm{d}\mathbf{W}_t`` with ``\mathbf{Q}_{\mathrm{user}} = \mathbf{\Sigma}_0\mathbf{\Sigma}_0^\top`` and want the action in *your* metric ``\mathbf{Q}_{\mathrm{user}}^{-1}``, multiply the returned action by ``\mathrm{tr}(\mathbf{Q}_{\mathrm{user}})/D``:
```math
\Phi_{FW}^{\mathbf{Q}_{\mathrm{user}}}[\phi] \;=\; \frac{\mathrm{tr}(\mathbf{Q}_{\mathrm{user}})}{D}\cdot
    \texttt{fw\_action(sys, path, time)}.
```
The factor is ``1`` whenever ``\mathbf{Q}_{\mathrm{user}}`` is isotropic (``c\,\mathbf{I}`` for any ``c>0``) or already trace-normalized (``\mathrm{tr}(\mathbf{Q}_{\mathrm{user}}) = D``), in which case `fw_action` returns ``\Phi_{FW}^{\mathbf{Q}_{\mathrm{user}}}`` directly. All examples and tests in the package fall in this case.

## Non-autonomous: `RateSystem`

```@docs; canonical=false
RateSystem
ForcingProfile
unforced_system
```

![Schematic explaining RateSystem construction](../assets/ratesystem_scheme.png)

### `RateSystem` API

```@docs; canonical=false
parameters
parameter
set_forcing_start!
set_forcing_duration!
set_forcing_scale!
```

## Converting between systems

The deterministic part of a [`CoupledSDEs`](@ref) can be extracted as a [`CoupledODEs`](@ref), making it compatible with functionality of `DynamicalSystems.jl`. In turn, a `CoupledODEs` can easily be extended into a `CoupledSDEs`.

```@docs; canonical=false
CoupledODEs(ds::CoupledSDEs; kwargs...)
CoupledSDEs(ds::CoupledODEs, p; kwargs...)
```

For example, the Lyapunov spectrum of a `CoupledSDEs` in the absence of noise, here exemplified by the FitzHugh-Nagumo model, can be computed by typing:

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