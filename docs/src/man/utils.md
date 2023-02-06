# Convenience functions and types

## Interface to `DynamicalSystems.jl`

!!! tip "Using the functionality of DynamicalSystems.jl"
    A [`StochSystem`](@ref) can easily be turned into a
    [`ContinuousDynamicalSystem`](https://juliadynamics.github.io/DynamicalSystems.jl/stable/ds/general/#Dynamical-System-Definition)
    of `DynamicalSystems.jl` using the function [`to_cds`](@ref).

This way, many of the methods in `DynamicalSystems.jl` can be used directly, even if we have
not written an analogue method that takes `StochSystem` as input. For example, the
Lyapunov spectrum of a `StochSystem`, here exemplified by the FitzHugh-Nagumo model, can be
computed by typing:

```julia
using CriticalTransitions, DynamicalSystems: lyapunovspectrum
sys = StochSystem(FitzHughNagumo, [1.,3.,1.,1.,1.,0.], 2)
ls = lyapunovspectrum(to_cds(sys), 10000)
```

```@docs
to_cds(sys::StochSystem; state=zeros(sys.dim))
tocds
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

### `sys.dim`-dimensional box

```@docs
intervals_to_box(bmin::Vector, bmax::Vector)
```

## Custom types

* `Parameters = Union{Vector{Float64}, Vector{Int64}, Nothing}`
* `CovMatrix = Union{Matrix, UniformScaling{Bool}, Diagonal{Bool, Vector{Bool}}}`
* `State = Union{Vector, SVector}`