module CriticalTransitions

# Base
using Statistics: Statistics, mean
using LinearAlgebra: LinearAlgebra, I, norm, dot, tr, det
using StaticArrays: StaticArrays, SVector
using SparseArrays: spdiagm
using DataStructures: CircularBuffer
using Random: Random

# Core
using SciMLBase: EnsembleThreads, DiscreteCallback, remake, terminate!
using StochasticDiffEq: StochasticDiffEq
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
using LinearSolve: LinearProblem, KLUFactorization, solve

# io and documentation
using Format: Format
using Printf: Printf
using DocStringExtensions: TYPEDSIGNATURES, TYPEDEF, TYPEDFIELDS, METHODLIST
using ProgressMeter: Progress, next!

# reexport
using Reexport: @reexport
@reexport using StaticArrays
@reexport using DynamicalSystemsBase

include("extension_functions.jl")
include("utils.jl")
include("sde_utils.jl")

include("trajectories/TransitionEnsemble.jl")
include("trajectories/simulation.jl")
include("trajectories/transition.jl")

include("transition_path_theory/langevin.jl")
include("transition_path_theory/committor.jl")
include("transition_path_theory/invariant_pdf.jl")
include("transition_path_theory/reactive_current.jl")
include("transition_path_theory/probability.jl")

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
export drift, div_drift, solver
export StateSpaceSet

export sgmam, SgmamSystem
export fw_action, om_action, action, geometric_action
export min_action_method, geometric_min_action_method, string_method
export MinimumActionPath

export deterministic_orbit
export transition, transitions

export basins, basinboundary, basboundary
export intervals_to_box

export distmesh2D, dellipse, ddiff
export get_ellipse, reparametrization
export find_boundary, huniform, dunion

export committor, invariant_pdf, reactive_current, probability_reactive, probability_last_A

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
