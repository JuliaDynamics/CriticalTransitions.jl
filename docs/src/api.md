# API

This page is the canonical reference for the public API, grouped by topic.
Each entry also appears, in context, on the relevant manual page; the canonical
link used by `@ref` points here.

## Rate tipping

### `RateSystem`

```@docs
RateSystem
ForcingProfile
unforced_system
parameters
parameter
set_forcing_start!
set_forcing_duration!
set_forcing_scale!
set_forcing_reverse!
```

### Rate tipping functions

```@docs
rate_track_return_tip
```

## Stochastic dynamical systems

```@docs
CoupledSDEs
CoupledSDEs(ds::CoupledODEs, p; kwargs...)
CoupledODEs(ds::CoupledSDEs; kwargs...)
FreidlinWentzellHamiltonian
```



### Accessors

```@docs
solver
drift
div_drift
covariance_matrix
diffusion_matrix
noise_process
noise_strength
```

## Simulation

```@docs
deterministic_orbit
```

## Sampling

```@docs
transition
transitions
CriticalTransitions.TransitionEnsemble
CriticalTransitions.TransitionStatistics
```

## Large deviations

### Action functionals

```@docs
fw_action
geometric_action
om_action
action
```

### Action minimizers

```@docs
minimize_action
minimize_geometric_action
string_method
CriticalTransitions.MinimumActionPath
```
#### Geometric minimal action methods

```@docs
GeometricGradient
AdaptiveGeometricGradient
MultipleShooting
```

### Quasipotential

```@docs
quasipotential
CriticalTransitions.QuasiPotential
CriticalTransitions.BackRef
```

## Generator and rates

### Generator and grid

```@docs
DiffusionGenerator
rate_matrix
m_matrix
fokker_planck_operator
CartesianGrid
CriticalTransitions.BoundaryCondition
Reflecting
Periodic
Absorbing
```

### Observables

```@docs
stationary_distribution
quasi_stationary_distribution
mean_first_passage_time
first_passage_variance
eigenmodes
propagate_density
```

### Eigensolvers

```@docs
DenseEigen
KrylovKitSolver
```

### Helpers

```@docs
CriticalTransitions.ncells
CriticalTransitions.cell_volume
CriticalTransitions.cell_center
CriticalTransitions.ball
CriticalTransitions.cuboid
CriticalTransitions.sublevel
CriticalTransitions.reshape_to_grid
```

## Reactive transitions (Transition path theory)

```@docs
ReactiveTransition
forward_committor
backward_committor
reactive_rate
reactive_density
reactive_current
reactive_current_reversible
reactive_current_irreversible
probability_reactive
probability_last_A
```

## Utilities

```@docs
CriticalTransitions.normalize_covariance!
```


