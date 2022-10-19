# Analyzing a system's stability properties

## Fixed points
```@docs
equilib(sys::StochSystem, state::State; kwargs...)
fixedpoints(sys::StochSystem, bmin, bmax)
```

## Basins of attraction
```@docs
basins(sys::StochSystem, A, B, C, H; kwargs)
basinboundary(X, Y, h; kwargs...)
```
