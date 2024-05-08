# Analyzing a system's stability properties

To use the following functionalities, you need to load `DynamicalSystems.jl`.

## Fixed points
```@docs
equilib(sys::CoupledSDEs, state::State; kwargs...)
fixedpoints
```

## Basins of attraction
```@docs
basins(sys::CoupledSDEs, A, B, C, H; kwargs)
basinboundary(X, Y, h; kwargs...)
basboundary(sys::CoupledSDEs, xrange::Vector, yrange::Vector, xspacing::Float64, attractors::Vector; kwargs...)
```

## Edge tracking
The edge tracking algorithm is a simple numerical method to find the *edge state* or
(possibly chaotic) saddle on the boundary between two basins of attraction. It is first
introduced by [Battelino et al. (1988)](https://doi.org/10.1016/0167-2789(88)90057-7) and further described by [Skufca et al. (2006)](https://doi.org/10.1103/PhysRevLett.96.174101).

```@docs
edgetracking(sys::CoupledSDEs, u1::State, u2::State, attractors::Vector; kwargs...)
bisect_to_edge(sys::CoupledSDEs, u1::State, u2::State, attractors::Vector; kwargs...)
```