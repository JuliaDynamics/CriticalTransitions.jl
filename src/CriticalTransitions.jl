module CriticalTransitions

# Use the README as the module docs
@doc let
    path = joinpath(dirname(@__DIR__), "README.md")
    include_dependency(path)
    read(path, String)
end CriticalTransitions

# Base
using Statistics: Statistics, mean
using SparseArrays: SparseArrays
using StaticArrays: StaticArrays, SVector, SMatrix, MVector
using LinearAlgebra: LinearAlgebra, norm, dot, tr, diag, eigen, normalize!, I
using Random: Random

# Core
using ForwardDiff: ForwardDiff
using SciMLBase: SciMLBase, EnsembleThreads, DiscreteCallback, remake, terminate!, isinplace
using OrdinaryDiffEqLowOrderRK: OrdinaryDiffEqLowOrderRK, Euler
using NonlinearSolveFirstOrder: NonlinearSolveFirstOrder
using ADTypes: ADTypes, AutoForwardDiff
using SparseMatrixColorings: SparseMatrixColorings, GreedyColoringAlgorithm
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
    set_parameter!,
    set_parameters!,
    initial_time,
    integrator,
    referenced_sciml_prob,
    covariance_matrix,
    diffusion_matrix,
    diffusion_function,
    current_time
using Attractors: Attractors

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
using DataStructures: DataStructures, MutableBinaryHeap, FasterForward

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
@reexport using Attractors

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
include("largedeviations/multiple_shooting.jl")
include("largedeviations/string_method.jl")

include("r_tipping/ForcingProfile.jl")
include("r_tipping/RateSystem.jl")
include("r_tipping/r_tipping_phase_diagrams.jl")

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
include("transition_path_theory/reactive_transition.jl")
include("largedeviations/quasipotential/state.jl")
include("largedeviations/quasipotential/stencil.jl")
include("largedeviations/quasipotential/lagrangian.jl")
include("largedeviations/quasipotential/update.jl")
include("largedeviations/quasipotential/sweep.jl")
include("largedeviations/quasipotential/quasipotential.jl")

include("../systems/CTLibrary.jl")
using .CTLibrary

# Core types
export CoupledSDEs, CoupledODEs, noise_process, covariance_matrix, diffusion_matrix
export drift, div_drift, solver, noise_strength
export StateSpaceSet

export FreidlinWentzellHamiltonian
export fw_action, om_action, action, geometric_action
export minimize_action, minimize_geometric_action, string_method
export MinimumActionPath
export GeometricGradient, AdaptiveGeometricGradient, MultipleShooting
export quasipotential

export deterministic_orbit
export transition, transitions

export ForcingProfile, RateSystem
export set_forcing_duration!, set_forcing_scale!, set_forcing_start!, set_forcing_reverse!
export unforced_system, parameters, parameter
export rate_track_return_tip, unforced_pcurve

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

# Transition path theory
export ReactiveTransition
export forward_committor, backward_committor
export reactive_rate, reactive_density, reactive_current
export reactive_current_reversible, reactive_current_irreversible
export probability_reactive, probability_last_A

end # module CriticalTransitions
