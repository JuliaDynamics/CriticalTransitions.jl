module CoupledSDEsBasin

using CriticalTransitions
using ChaosTools
using Attractors:
    Attractors, AttractorsViaProximity, AttractorMapper, edgetracking, bisect_to_edge
using ProgressMeter
using DynamicalSystemsBase: ParallelDynamicalSystem
using DocStringExtensions
using LinearAlgebra

include("basin/planeofbox.jl")
include("basin/basinsofattraction.jl")
include("basin/edgetrack.jl")
include("basin/basinboundary.jl")

export basins, basboundary, basinboundary

end # module CoupledSDEsBasin
