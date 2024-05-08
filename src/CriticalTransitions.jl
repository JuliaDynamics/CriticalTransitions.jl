module CriticalTransitions

using Reexport
@reexport using DynamicalSystemsBase
@reexport using StaticArrays
@reexport using OrdinaryDiffEq
@reexport using StochasticDiffEq
@reexport using DiffEqNoiseProcess
@reexport using LinearAlgebra
using Format, Dates, JLD2, HDF5, ProgressBars, ProgressMeter, DocStringExtensions
using Attractors
using ChaosTools
using IntervalRootFinding
using ForwardDiff
using Symbolics
using Optim, Dierckx
using Printf, DrWatson, Dates, Statistics

include("utils.jl")
include("CoupledSDEs.jl")
include("io.jl")
include("trajectories/simulation.jl")
include("trajectories/transition.jl")
include("noiseprocesses/stochprocess.jl")
include("largedeviations/action.jl")
include("largedeviations/min_action_method.jl")
include("largedeviations/geometric_min_action_method.jl")

include("../systems/fitzhughnagumo.jl")
include("../systems/truscottbrindley_mod.jl")
include("../systems/truscottbrindley_orig.jl")
include("../systems/truscottbrindley_orig1.jl")
include("../systems/rooth.jl")
include("../systems/stommel.jl")
include("../systems/rivals.jl")

# include("../dev/fhn_pathspace_sampling.jl")
# include("../dev/symbolic_langevinmcmc.jl")
# include("../dev/residence_times.jl")
# include("../dev/edgetrack_ct.jl")
# include("../dev/flexibletransitions.jl")
# include("../dev/RateSys1.jl")

# Core types
export CoupledSDEs, diag_noise_funtion
export noise_strength

# Methods
export equilib, fixedpoints, basins, basinboundary, basboundary
export simulate, relax
export transition, transitions
export fw_integrand, fw_action, om_action, action, geometric_action
export min_action_method, geometric_min_action_method
# export edgetracking, bisect_to_edge, attractor_mapper, bisect_to_edge2
export make_jld2, make_h5, intervals_to_box

# Systems
export fitzhugh_nagumo, fitzhugh_nagumo!, fhn_ϵσ, fhn_ϵσ_backward
export modifiedtruscottbrindley, modifiedtruscottbrindley!, modtb_αξσ, modtb_αξσ1, modtb_αξσ_backward
export rampedmodifiedtruscottbrindley, modifiedtruscottbrindley!, rmodtb_ξvTtrTraσ
export originaltruscottbrindley, originaltruscottbrindley!, origtb_rσ
export rampedoriginaltruscottbrindley, rampedoriginaltruscottbrindley!, rorigtb_vTtrTraσ
export originaltruscottbrindley1, originaltruscottbrindley1!, origtb1_rσ
export rampedoriginaltruscottbrindley1, rampedoriginaltruscottbrindley1!, rorigtb1_vTtrTraσ
export rivals!, rivals, rivals_ϵσ
export rooth_smooth, stommel, cessi

end # module CriticalTransitions
