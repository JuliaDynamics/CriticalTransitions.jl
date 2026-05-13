module CriticalTransitions

# Use the README as the module docs
@doc let
    path = joinpath(dirname(@__DIR__), "README.md")
    include_dependency(path)
    read(path, String)
end CriticalTransitions

# Base
using Statistics: Statistics, mean
using LinearAlgebra:
    LinearAlgebra, norm, dot, tr, det, diag, eigen, normalize!, I
using StaticArrays: StaticArrays, SVector
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
    integrator,
    referrenced_sciml_prob

using ConstructionBase: ConstructionBase
using StateSpaceSets: StateSpaceSets, dimension, StateSpaceSet
using StochasticDiffEq: StochasticDiffEq

using SparseArrays:
    SparseArrays, sparse, SparseMatrixCSC, rowvals, nonzeros, nzrange

using Interpolations: linear_interpolation
using Optimization: Optimization
using OptimizationOptimisers: Optimisers
using LinearSolve:
    LinearProblem, LUFactorization, UMFPACKFactorization, init, solve, solve!
using KrylovKit: KrylovKit
using ExponentialUtilities: expv_timestep

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

# Diffusion operator (general SDE machinery: discrete generator + analyses)
include("diffusion_operator/utils.jl")
include("diffusion_operator/cartesian_grid.jl")
include("diffusion_operator/diffusion_generator.jl")
include("diffusion_operator/spectral.jl")
include("diffusion_operator/generator_analyses.jl")
include("diffusion_operator/propagator.jl")
include("diffusion_operator/grid_helpers.jl")

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

# Diffusion operator
export CartesianGrid, DiffusionGenerator
export Reflecting, Periodic, Absorbing
export rate_matrix, m_matrix, fokker_planck_operator
export stationary_distribution
export quasi_stationary_distribution
export mean_first_passage_time, first_passage_variance
export eigenmodes
export propagate_density
export DenseEigen, KrylovKitSolver

end # module CriticalTransitions
