module CriticalTransitions

# Use the README as the module docs
@doc let
    path = joinpath(dirname(@__DIR__), "README.md")
    include_dependency(path)
    read(path, String)
end CriticalTransitions

# Base
using Statistics: Statistics, mean
using LinearAlgebra: LinearAlgebra, norm, dot, tr, det
using StaticArrays: StaticArrays, SVector
using DataStructures: CircularBuffer
using Random: Random

# Core
using SciMLBase: SciMLBase, EnsembleThreads, DiscreteCallback, remake, terminate!, isinplace
using OrdinaryDiffEqLowOrderRK: OrdinaryDiffEqLowOrderRK, Euler
using DynamicalSystemsBase:
    DynamicalSystemsBase,
    CoupledSDEs,
    CoupledODEs,
    dynamic_rule,
    initial_state,
    current_state,
    set_state!,
    trajectory,
    jacobian,
    ContinuousTimeDynamicalSystem,
    initial_parameters,
    current_parameter,
    current_parameters,
    set_parameters!,
    initial_time,
    integrator

using ConstructionBase: ConstructionBase
using StateSpaceSets: StateSpaceSets, dimension, StateSpaceSet
using StochasticDiffEq: StochasticDiffEq

using Interpolations: linear_interpolation
using Optimization: Optimization
using OptimizationOptimisers: Optimisers
using LinearSolve: LinearProblem, LUFactorization, init, solve, solve!

# io and documentation
using Format: Format
using Printf: Printf
using DocStringExtensions:
    TYPEDSIGNATURES, TYPEDEF, TYPEDFIELDS, METHODLIST, DocStringExtensions
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

include("largedeviations/utils.jl")
include("largedeviations/action.jl")
include("largedeviations/methods.jl")
include("largedeviations/MinimumActionPath.jl")
include("largedeviations/minimize_action.jl")
include("largedeviations/minimize_geometric_action.jl")

include("largedeviations/sgMAM.jl")
include("largedeviations/string_method.jl")

include("r_tipping/RateSystem.jl")

# Experimental features
include("experimental/transition_path_theory/TransitionPathMesh.jl")
include("experimental/transition_path_theory/langevin.jl")
include("experimental/transition_path_theory/committor.jl")
include("experimental/transition_path_theory/invariant_pdf.jl")
include("experimental/transition_path_theory/reactive_current.jl")
include("experimental/transition_path_theory/probability.jl")

include("../systems/CTLibrary.jl")
using .CTLibrary

# Core types
export CoupledSDEs, CoupledODEs, noise_process, covariance_matrix, diffusion_matrix
export drift, div_drift, solver
export StateSpaceSet

export minimize_simple_geometric_action, ExtendedPhaseSpace
export fw_action, om_action, action, geometric_action
export minimize_action, action_minimizer, minimize_geometric_action, string_method
export MinimumActionPath
export GeometricGradient

export deterministic_orbit
export transition, transitions

export ForcingProfile, RateSystem
export set_forcing_duration!, set_forcing_scale!, set_forcing_start!
export frozen_system, parameters

end # module CriticalTransitions
