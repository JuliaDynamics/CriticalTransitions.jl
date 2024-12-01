module CriticalTransitions

# Base
using Statistics: Statistics, mean
using LinearAlgebra: LinearAlgebra, I, norm, dot, tr
using StaticArrays: StaticArrays, SVector

# Core
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
using DynamicalSystemsBase:
    DynamicalSystemsBase,
    CoupledSDEs,
    CoupledODEs,
    dynamic_rule,
    current_state,
    set_state!,
    trajectory

using ForwardDiff: ForwardDiff
using IntervalArithmetic: IntervalArithmetic, interval
using Interpolations: linear_interpolation
using Optim: Optim, LBFGS
using Symbolics: Symbolics

# io and documentation
using Format: Format
using Dates: Dates
using Printf: Printf
using Markdown: Markdown
using DocStringExtensions: TYPEDSIGNATURES
using ProgressBars: ProgressBars, tqdm
using ProgressMeter: ProgressMeter

# reexport
using Reexport: @reexport
@reexport using StaticArrays
@reexport using StochasticDiffEq
@reexport using DiffEqNoiseProcess

include("extension_functions.jl")
include("utils.jl")
include("system_utils.jl")
include("trajectories/simulation.jl")
include("trajectories/transition.jl")
include("trajectories/equib.jl")
include("noiseprocesses/stochprocess.jl")
include("largedeviations/action.jl")
include("largedeviations/min_action_method.jl")
include("largedeviations/geometric_min_action_method.jl")

include("largedeviations/sgMAM.jl")
using .Sgmam: sgmam, SgmamSystem

include("../systems/CTLibrary.jl")
using .CTLibrary

# Core types
export CoupledSDEs, CoupledODEs, noise_process, covariance_matrix, diffusion_matrix
export dynamic_rule, current_state, set_state!, trajectory

export sgmam, SgmamSystem

# Methods
export drift, div_drift
export equilib, deterministic_orbit
export transition, transitions
export basins, basinboundary, basboundary
export fw_action, om_action, action, geometric_action
export min_action_method, geometric_min_action_method
export intervals_to_box
export covariance_matrix, diffusion_matrix
# export edgetracking, bisect_to_edge, AttractorsViaProximity
# export fixedpoints
# ^ extension tests needed

# Error hint for extensions stubs
function __init__()
    Base.Experimental.register_error_hint(_basin_error_hinter(basins), MethodError)
    Base.Experimental.register_error_hint(_basin_error_hinter(basboundary), MethodError)
    Base.Experimental.register_error_hint(_basin_error_hinter(basinboundary), MethodError)
    return nothing
end

end # module CriticalTransitions
