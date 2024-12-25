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
using Interpolations: linear_interpolation
using Optim: Optim, LBFGS
using Optimization
using OptimizationOptimisers: Optimisers
using Symbolics: Symbolics

# io and documentation
using Format: Format
using Printf: Printf
using DocStringExtensions: TYPEDSIGNATURES
using ProgressMeter: Progress, next!

# reexport
using Reexport: @reexport
@reexport using StaticArrays
@reexport using StochasticDiffEq
@reexport using DiffEqNoiseProcess

include("extension_functions.jl")
include("utils.jl")
include("sde_utils.jl")
include("trajectories/simulation.jl")
include("trajectories/transition.jl")
include("trajectories/equib.jl")
include("noiseprocesses/stochprocess.jl")
include("largedeviations/action.jl")
include("largedeviations/MaximumLikelihoodPath.jl")
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
export intervals_to_box
export fw_action, om_action, action, geometric_action
export min_action_method, geometric_min_action_method
export covariance_matrix, diffusion_matrix
# export edgetracking, bisect_to_edge, AttractorsViaProximity

# Error hint for extensions stubs
function __init__()
    Base.Experimental.register_error_hint(_basin_error_hinter(basins), MethodError)
    Base.Experimental.register_error_hint(_basin_error_hinter(basboundary), MethodError)
    Base.Experimental.register_error_hint(_basin_error_hinter(basinboundary), MethodError)
    Base.Experimental.register_error_hint(
        _basin_error_hinter(intervals_to_box), MethodError
    )
    return nothing
end

end # module CriticalTransitions
