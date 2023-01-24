module CriticalTransitions

using Formatting, Dates, JLD2, HDF5, ProgressBars, ProgressMeter, DocStringExtensions
using DynamicalSystems, IntervalRootFinding
using OrdinaryDiffEq, StochasticDiffEq, DiffEqNoiseProcess
using LinearAlgebra, StaticArrays, ForwardDiff

include("StochSystem.jl")
include("utils.jl")

include("io/io.jl")
include("noiseprocesses/gaussian.jl")
include("systemanalysis/stability.jl")
include("systemanalysis/basinsofattraction.jl")
include("systemanalysis/basinboundary.jl")
include("trajectories/simulation.jl")
include("trajectories/transition.jl")
include("pathproperties/action.jl")

include("../systems/fitzhughnagumo.jl")
include("../systems/modifiedtruscottbrindley.jl")

include("../dev/fhn_pathspace_sampling.jl")

export StochSystem, State
export equilib, fixedpoints, basins, basinboundary
export simulate, relax, transition, transitions
export tocds, make_jld2, make_h5, sys_string, sys_info, intervals_to_box
export is_iip
export fw_integrand, fw_action

export FitzHughNagumo, FitzHughNagumo!
export modifiedtruscottbrindley, modifiedtruscottbrindley!
export modtb_αξσ
export idfunc, idfunc!
export additive_idx, additive_idx!
export multiplicative_idx, multiplicative_idx!
export gauss

export EM, I # functions inherited from dependencies

export FitzHughNagumoSPDE, fhn_pathspace_sampling

end # module CriticalTransitions