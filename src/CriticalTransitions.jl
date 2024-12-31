module CriticalTransitions

# Base
using Statistics: Statistics, mean
using LinearAlgebra: LinearAlgebra, I, norm, dot, tr
using StaticArrays: StaticArrays, SVector
using SparseArrays: spdiagm
using DataStructures: CircularBuffer

# Core
using DiffEqNoiseProcess: DiffEqNoiseProcess
using OrdinaryDiffEq: OrdinaryDiffEq, EnsembleThreads
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
    trajectory,
    jacobian,
    StateSpaceSet

using Interpolations: linear_interpolation
using Optimization
using OptimizationOptimisers: Optimisers
using Symbolics: Symbolics
using LinearSolve: LinearProblem, KLUFactorization, solve

# io and documentation
using Format: Format
using Printf: Printf
using DocStringExtensions: TYPEDSIGNATURES, TYPEDEF, TYPEDFIELDS, METHODLIST
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

include("largedeviations/utils.jl")
include("largedeviations/action.jl")
include("largedeviations/MinimumActionPath.jl")
include("largedeviations/min_action_method.jl")
include("largedeviations/geometric_min_action_method.jl")

include("largedeviations/sgMAM.jl")
include("largedeviations/string_method.jl")

include("../systems/CTLibrary.jl")
using .CTLibrary

# Core types
export CoupledSDEs, CoupledODEs, noise_process, covariance_matrix, diffusion_matrix
export dynamic_rule, current_state, set_state!, trajectory
export drift, div_drift, solver
export StateSpaceSet

export sgmam, SgmamSystem
export fw_action, om_action, action, geometric_action
export min_action_method, geometric_min_action_method, string_method
export MinimumActionPath

export deterministic_orbit
export transition, transitions
export TransitionEnsemble

export basins, basinboundary, basboundary
export intervals_to_box

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
