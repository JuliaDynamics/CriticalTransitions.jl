# Convenience functions and types

## Interface to `DynamicalSystems.jl`

!!! tip "Using the functionality of DynamicalSystems.jl"
    A [`CoupledSDEs`](@ref) can easily be turned into a 
    [`CoupledODEs`](https://juliadynamics.github.io/DynamicalSystems.jl/dev/tutorial/#DynamicalSystemsBase.CoupledODEs)
    instance of `DynamicalSystems.jl` using the function [`CoupledODEs`](@ref). Vice vera,
    a `CoupledODEs` system can be converted into a `CoupledSDEs` via `CoupledSDEs(ds::CoupledODEs, g)` with g the noise function.

This way, many of the methods in `DynamicalSystems.jl` can be used directly, even if we have
not written an analogue method that takes `CoupledSDEs` as input. For example, the
Lyapunov spectrum of a `CoupledSDEs`, here exemplified by the FitzHugh-Nagumo model, can be
computed by typing:

```julia
using CriticalTransitions, DynamicalSystems: lyapunovspectrum
sys = CoupledSDEs(fitzhugh_nagumo, diag_noise_function(σ), zeros(2), [1.,3.,1.,1.,1.,0.])
ls = lyapunovspectrum(CoupledODEs(sys), 10000)
```

```@docs
CoupledODEs(sys::CoupledSDEs; diffeq, t0=0.0)
```

## `CoupledSDEs` utility functions

```@docs
CriticalTransitions.drift
```

## Saving data

```@docs
make_jld2(text::String, relpath::String="")
make_h5(text::String, relpath::String="")
```

### `length(sys.u)`-dimensional box

```@docs
intervals_to_box(bmin::Vector, bmax::Vector)
```
