module CriticalTransitions

using Formatting, Dates, JLD2, HDF5, ProgressBars, ProgressMeter, DocStringExtensions
using DynamicalSystems, IntervalRootFinding
using OrdinaryDiffEq, StochasticDiffEq, DiffEqNoiseProcess
using LinearAlgebra, StaticArrays, ForwardDiff
using Symbolics
using Optim, Dierckx
using Printf, DrWatson, Dates, Statistics

include("StochSystem.jl")
include("utils.jl")

include("io/io.jl")
include("noiseprocesses/gaussian.jl")
include("systemanalysis/stability.jl")
include("systemanalysis/basinsofattraction.jl")
include("systemanalysis/basinboundary.jl")
include("trajectories/simulation.jl")
include("trajectories/transition.jl")
include("largedeviations/action.jl")
include("largedeviations/mam.jl")
include("largedeviations/gmam.jl")

include("../systems/fitzhughnagumo.jl")
include("../systems/modifiedtruscottbrindley.jl")
include("../systems/thc_box/rooth_model.jl")
include("../systems/rivals.jl")

include("../dev/fhn_pathspace_sampling.jl")
include("../dev/symbolic_langevinmcmc.jl")
include("../dev/residence_times.jl")
include("../dev/edgetrack_ct.jl")
include("../dev/flexibletransitions.jl")

export StochSystem, State
export equilib, fixedpoints, basins, basinboundary
export saddles_idx, repellers_idx, attractors_idx
export simulate, relax, transition, transitions
export tocds, to_cds, drift
export make_jld2, make_h5, sys_string, sys_info, intervals_to_box
export is_iip
export fw_integrand, fw_action, om_action, geometric_action
export mam, gmam

export FitzHughNagumo, FitzHughNagumo!, fhn_ϵσ
export modifiedtruscottbrindley, modifiedtruscottbrindley!, modtb_αξσ
export rivals!, rivals, rivals_ϵσ 
export rooth_smooth
export transition2, transitions2
export residence_time2, residence_times2

export idfunc, idfunc!
export additive_idx, additive_idx!
export multiplicative_idx, multiplicative_idx!
export anorm, subnorm
export gauss

export EM, I # functions inherited from dependencies

export FitzHughNagumoSPDE, fhn_pathspace_sampling
export langevinmcmc_spde, symbolise_spde, langevinmcmc
export jacobian
export residence_time, residence_times, ResTimes, temporal, runandsavetimes, get_res_times
export edgetracking, bisect_to_edge, attractor_mapper

end # module CriticalTransitions