# Analyzing a system's stability properties

## Fixed points
```@docs
equilib(sys::StochSystem, state::State; kwargs...)
fixedpoints
```

## Basins of attraction
```@docs
basins(sys::StochSystem, A, B, C, H; kwargs)
basinboundary(X, Y, h; kwargs...)
```

## Edge tracking
The edge tracking algorithm is a simple numerical method to find the *edge state* or
(possibly chaotic) saddle on the boundary between two basins of attraction. It is first
described by [Skufca et al. (2006)](https://journals.aps.org/prl/pdf/10.1103/PhysRevLett.96.174101?casa_token=RUn26KnFdNEAAAAA%3AoXsTlmEWVMkEYbOtR-j2PH2vYOOPOy1a2R_37ncnf4gsiHp6GR66M-IBpzXocLoMQC_oHhk8MIFRa_8).

```@docs
edgetracking(sys::StochSystem, u1::State, u2::State, attractors::Vector; kwargs...)
bisect_to_edge(sys::StochSystem, u1::State, u2::State, attractors::Vector; kwargs...)
```