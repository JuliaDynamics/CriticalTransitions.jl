module CriticalTransitions

using Reexport
@reexport using DynamicalSystemsBase
@reexport using StaticArrays
@reexport using OrdinaryDiffEq
@reexport using StochasticDiffEq
@reexport using DiffEqNoiseProcess
@reexport using LinearAlgebra
using Format, Dates, JLD2, HDF5, ProgressBars, ProgressMeter, DocStringExtensions
using IntervalRootFinding
using ForwardDiff
using Symbolics
using Optim, Dierckx
using Printf, DrWatson, Dates, Statistics

# Define custom types
# Parameters = Union{Vector{Any}, Nothing};
# CovMatrix = Union{Matrix, UniformScaling{Bool}, Diagonal{Bool, Vector{Bool}}};
# State = Union{Vector, SVector}

include("extention_functions.jl")
include("utils.jl")
include("CoupledSDEs.jl")
include("io.jl")
include("trajectories/simulation.jl")
include("trajectories/transition.jl")
include("trajectories/equib.jl")
include("noiseprocesses/stochprocess.jl")
include("largedeviations/action.jl")
include("largedeviations/min_action_method.jl")
include("largedeviations/geometric_min_action_method.jl")

include("../systems/CTLibrary.jl")

# Core types
export CoupledSDEs, diagonal_noise!, diagonal_noise, add_noise_strength

# Methods
export equilib, fixedpoints, basins, basinboundary, basboundary
export simulate, relax
export transition, transitions
export fw_integrand, fw_action, om_action, action, geometric_action
export min_action_method, geometric_min_action_method
export basins, basinboundary
export edgetracking, bisect_to_edge, attractor_mapper
export make_jld2, make_h5, intervals_to_box

end # module CriticalTransitions
