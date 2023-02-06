# Convenience functions and types

## Functions

### `StochSystem` utility functions

```@docs
tocds(sys::StochSystem; state=zeros(sys.dim))
sys_info(sys::StochSystem)
sys_string(sys::StochSystem; verbose=true)
```

### Noise functions for `sys.g`

```@docs
idfunc(u, p, t)
idfunc!(du, u, p, t)
```

### Vector norms
```@docs
anorm(vec, A; square=false)
subnorm(vec; directions=[1,...,N])
```

### `sys.dim`-dimensional box

```@docs
intervals_to_box(bmin::Vector, bmax::Vector)
```

### Saving data

```@docs
make_jld2(text::String, relpath::String="")
make_h5(text::String, relpath::String="")
```

### Miscellaneous

```@docs
is_iip(f::Function)
```

## Types

* `Parameters = Union{Vector{Float64}, Vector{Int64}, Nothing}`
* `CovMatrix = Union{Matrix, UniformScaling{Bool}, Diagonal{Bool, Vector{Bool}}}`
* `State = Union{Vector, SVector}`