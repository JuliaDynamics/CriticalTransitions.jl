module CriticalTransitions

using Formatting, Dates, JLD2, HDF5, ProgressBars, ProgressMeter, DocStringExtensions
using DynamicalSystems, IntervalRootFinding
using OrdinaryDiffEq, StochasticDiffEq, DiffEqNoiseProcess
using LinearAlgebra, StaticArrays, ForwardDiff
using Symbolics
using Optim
using Printf, DrWatson, Dates

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
include("actionminimizers/mam.jl")

include("../systems/fitzhughnagumo.jl")
include("../systems/modifiedtruscottbrindley.jl")
include("../systems/thc_box/rooth_model.jl")

include("../dev/fhn_pathspace_sampling.jl")
include("../dev/symbolic_langevinmcmc.jl")
include("../dev/residence_times.jl")

export StochSystem, State
export equilib, fixedpoints, basins, basinboundary
export simulate, relax, transition, transitions
export tocds, make_jld2, make_h5, sys_string, sys_info, intervals_to_box
export is_iip
export fw_integrand, fw_action
export mam

export FitzHughNagumo, FitzHughNagumo!, fhn_ϵσ
export modifiedtruscottbrindley, modifiedtruscottbrindley!, modtb_αξσ
export rooth_smooth

export idfunc, idfunc!
export additive_idx, additive_idx!
export multiplicative_idx, multiplicative_idx!
export gauss

export EM, I # functions inherited from dependencies

export FitzHughNagumoSPDE, fhn_pathspace_sampling
export langevinmcmc_spde, symbolise_spde, langevinmcmc
export residence_time, residence_times, ResTimes, temporal, runandsavetimes

end # module CriticalTransitions