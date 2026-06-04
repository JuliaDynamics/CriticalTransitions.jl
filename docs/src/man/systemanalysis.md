# Analyzing a system's stability properties

To use the following functionalities, you need to load `ChoasTools.jl` and `Attractors.jl` (which are included as dependencies in CriticalTransitions.jl).

## Fixed points

The fixed points of the deterministic part of a [`CoupledSDEs`](@ref) can be found with [`ChaosTools.fixedpoints`](@extref); see the ChaosTools documentation for details.

## Edge tracking
The edge tracking algorithm is a simple numerical method to find the *edge state* or (possibly chaotic) saddle on the boundary between two basins of attraction. It is first introduced by [Battelino1988](@citet) and further described by [Skufca2006, Mehling2024](@citet).

See [`Attractors.edgetracking`](@extref), [`Attractors.bisect_to_edge`](@extref), and [`Attractors.EdgeTrackingResults`](@extref) in the Attractors documentation.