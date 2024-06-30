module CoupledSDEsBaisin

using CriticalTransitions
using ChaosTools
using Attractors: Attractors, AttractorsViaProximity, AttractorMapper, edgetracking, bisect_to_edge
using ProgressMeter
using ProgressBars
using DynamicalSystemsBase: ParallelDynamicalSystem
using DocStringExtensions

include("basin/planeofbox.jl")
include("basin/basinsofattraction.jl")
include("basin/edgetrack.jl")
include("basin/basinboundary.jl")

export basins, basboundary, basinboundary

end # module CoupledSDEsBaisin
