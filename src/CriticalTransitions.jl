module CriticalTransitions

# Base
using Statistics: Statistics, mean
using LinearAlgebra: LinearAlgebra, Diagonal, I, norm, tr, dot
using StaticArrays: StaticArrays, SVector

# core
using DynamicalSystemsBase:
    DynamicalSystemsBase,
    ContinuousTimeDynamicalSystem,
    StateSpaceSets,
    dimension,
    dynamic_rule
using DiffEqNoiseProcess: DiffEqNoiseProcess
using OrdinaryDiffEq: OrdinaryDiffEq, Tsit5
using StochasticDiffEq:
    StochasticDiffEq,
    DiscreteCallback,
    ODEProblem,
    SDEFunction,
    SOSRA,
    remake,
    solve,
    step!,
    terminate!,
    u_modified!

using ForwardDiff: ForwardDiff
using IntervalArithmetic: IntervalArithmetic, interval
using Interpolations: LinearInterpolation
using Optim: Optim, LBFGS
using Symbolics: Symbolics

# io and documentation
using Format: Format
using Dates: Dates
using Printf: Printf
using Markdown: Markdown
using DocStringExtensions: TYPEDSIGNATURES
using HDF5: HDF5, h5open, push!
using JLD2: JLD2, jldopen
using ProgressBars: ProgressBars, tqdm
using ProgressMeter: ProgressMeter

# reexport
using Reexport: @reexport
@reexport using StaticArrays
@reexport using StochasticDiffEq
@reexport using DiffEqNoiseProcess

include("extention_functions.jl")
include("utils.jl")
include("CoupledSDEs.jl")
include("CoupledSDEs_utils.jl")
include("io.jl")
include("trajectories/simulation.jl")
include("trajectories/transition.jl")
include("trajectories/equib.jl")
include("noiseprocesses/stochprocess.jl")
include("largedeviations/action.jl")
include("largedeviations/min_action_method.jl")
include("largedeviations/geometric_min_action_method.jl")

include("../systems/CTLibrary.jl")
using .CTLibrary

# Core types
export CoupledSDEs,
    idfunc!,
    idfunc,
    add_noise_strength,
    noise_process,
    covariance_matrix,
    noise_strength,
    CoupledODEs

# Methods
export equilib, basins, basinboundary, basboundary
export simulate, relax
export transition, transitions
export fw_action, om_action, action, geometric_action
export min_action_method, geometric_min_action_method
export make_jld2, make_h5, intervals_to_box
# export basins, basinboundary
# export edgetracking, bisect_to_edge, AttractorsViaProximity
# export fixedpoints
# ^ extention tests needed
end # module CriticalTransitions
