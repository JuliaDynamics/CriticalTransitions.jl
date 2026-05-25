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
using SparseArrays: SparseArrays
using StaticArrays: StaticArrays, SVector
using Random: Random

# Core
using ForwardDiff: ForwardDiff
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
    referenced_sciml_prob,
    covariance_matrix,
    diffusion_matrix,
    diffusion_function

using ConstructionBase: ConstructionBase
using StateSpaceSets: StateSpaceSets, dimension, StateSpaceSet
using StochasticDiffEq: StochasticDiffEq

using SparseArrays:
    SparseArrays, sparse, SparseMatrixCSC, rowvals, nonzeros, nzrange

using OptimizationBase: OptimizationBase
using OptimizationOptimisers: Optimisers
using KrylovKit: KrylovKit
using ExponentialUtilities: expv_timestep
using FastInterpolations: linear_interp!
using LinearSolve: LinearSolve, LinearProblem, LUFactorization, UMFPACKFactorization, init, solve, solve!

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

include("utils.jl")
include("sde_utils.jl")

include("trajectories/TransitionEnsemble.jl")
include("trajectories/simulation.jl")
include("trajectories/transition.jl")

include("largedeviations/utils.jl")
include("largedeviations/hamiltonian.jl")
include("largedeviations/action.jl")
include("largedeviations/methods.jl")
include("largedeviations/MinimumActionPath.jl")
include("largedeviations/minimize_action.jl")
include("largedeviations/sgmam_kernels.jl")
include("largedeviations/sgmam.jl")
include("largedeviations/minimize_geometric_action.jl")
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

# Transition path theory on the discrete diffusion generator
include("transition_path_theory/committor.jl")

include("../systems/CTLibrary.jl")
using .CTLibrary

# Core types
export CoupledSDEs, CoupledODEs, noise_process, covariance_matrix, diffusion_matrix
export drift, div_drift, solver, noise_strength
export StateSpaceSet

export FreidlinWentzellHamiltonian
export fw_action, om_action, action, geometric_action
export minimize_action, action_minimizer, minimize_geometric_action, string_method
export MinimumActionPath
export GeometricGradient, AdaptiveGeometricGradient

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
