# Convenience functions and types

## Interface to `DynamicalSystems.jl`

!!! tip "Using the functionality of DynamicalSystems.jl"
    A [`StochSystem`](@ref) can easily be turned into a 
    [`CoupledODEs`](https://juliadynamics.github.io/DynamicalSystems.jl/stable/ds/general/#Dynamical-System-Definition)
    instance of `DynamicalSystems.jl` using the function [`CoupledODEs`](@ref). Vice vera,
    a `CoupledODEs` system can be converted into a `StochSystem` via `StochSystem(ds::CoupledODEs)`.

This way, many of the methods in `DynamicalSystems.jl` can be used directly, even if we have
not written an analogue method that takes `StochSystem` as input. For example, the
Lyapunov spectrum of a `StochSystem`, here exemplified by the FitzHugh-Nagumo model, can be
computed by typing:

```julia
using CriticalTransitions, DynamicalSystems: lyapunovspectrum
sys = StochSystem(fitzhugh_nagumo, [1.,3.,1.,1.,1.,0.], zeros(2))
ls = lyapunovspectrum(CoupledODEs(sys), 10000)
```

```@docs
CoupledODEs(sys::StochSystem; diffeq, t0=0.0)
StochSystem(ds::CoupledODEs, σ, g, pg, Σ, process)
attractor_mapper(sys::StochSystem, attractors, eps=0.01; kwargs...)
```

## `StochSystem` utility functions

```@docs
drift(sys::StochSystem, x::State)
sys_info(sys::StochSystem)
sys_string(sys::StochSystem; verbose=true)
```

## Noise functions for `sys.g`

```@docs
idfunc(u, p, t)
idfunc!(du, u, p, t)
```

## Saving data

```@docs
make_jld2(text::String, relpath::String="")
make_h5(text::String, relpath::String="")
```

## Miscellaneous

```@docs
is_iip(f::Function)
```

### Generalized vector norms
```@docs
anorm(vec, A; square=false)
subnorm(vec; directions=[1,...,N])
```

### `length(sys.u)`-dimensional box

```@docs
intervals_to_box(bmin::Vector, bmax::Vector)
```

## Custom types

* `Parameters = Union{Vector{Float64}, Vector{Int64}, Nothing}`
* `CovMatrix = Union{Matrix, UniformScaling{Bool}, Diagonal{Bool, Vector{Bool}}}`
* `State = Union{Vector, SVector}`