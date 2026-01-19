# Changelog for `CriticalTransitions.jl`

## v0.7.0
New `RateSystem` type

Introduced the following new types:
- `RateSystem`
- `ForcingProfile`
- `RateSystemSpecs` (behind the scenes data structure; not exported)

... along with several interface methods.

Several small other changes and revised docs. Full changelog [here](https://github.com/JuliaDynamics/CriticalTransitions.jl/compare/v0.6.0...v0.7.0).

## v0.6.0
New features and examples

- Several breaking changes due to renaming of functions and types.
- Functionality on transition path theory has been revised. Associated functions are not automatically exported but can be accessed by loading them explicitly.

Full changelog [here](https://github.com/JuliaDynamics/CriticalTransitions.jl/compare/v0.5.0...v0.6.0).

## v0.5.0
Package clean-up

This release brings a streamlined source code, with unused and deprecated code removed, dependencies reduced, and a unified style of displaying progress meters.

#### Changed functions
- `intervals_to_box` is now part of the `ChaosToolsExt` extension and is loaded only when `using ChaosTools`
- The keyword argument `showprogress` changed to `show_progress` in several functions for consistency
- In `min_action_method`, the keyword argument `Stol` is now called `action_tol`

#### Deprecations
- Functions `make_h5` and `make_jld2` have been removed
- Dependencies `HDF5`, `JLD2`, and `ProgressBars` have been removed
- Dependency `IntervalArithmetic` is now only a dependency of some extensions

Full changelog [here](https://github.com/JuliaDynamics/CriticalTransitions.jl/compare/v0.4.0...v0.5.0).

## v0.4.0
New `CoupledSDEs` design

This release upgrades CriticalTransitions.jl to be compatible with the re-design of `CoupledSDEs`, which has now been integrated in [`DynamicalSystemsBase.jl v3.11`](https://juliadynamics.github.io/DynamicalSystemsBase.jl/stable/CoupledSDEs/).

Since we have updated the syntax of defining a `CoupledSDEs` system, this is a breaking change.

#### Changed functions
- `CoupledSDEs` is now constructed with different args and kwargs
- `fw_action`, `geometric_action` and `om_action` now normalize the covariance matrix when computing the action functional
- `noise_strength` is not a function anymore but a kwarg for `CoupledSDEs` and other functions
- `om_action` now requires `noise_strength` as input argument
- `relax` was renamed to `deterministic_orbit`
- `trajectory` has been added to replace `simulate`

#### Deprecations
- `add_noise_strength`, `idfunc` and `idfunc!` are no longer needed and have been removed
- the function `noise_strength(::CoupledSDEs)` has been removed
- `simulate` has been removed

Full changelog [here](https://github.com/JuliaDynamics/CriticalTransitions.jl/compare/v0.3.0...v0.4.0).

## v0.3.0
Major overhaul introducing `CoupledSDEs`

This release replaces the `StochSystem` struct with the new `CoupledSDEs` struct to define a stochastic dynamical system. To see how the new version works, check out the [documentation](https://juliadynamics.github.io/CriticalTransitions.jl/dev/).

The update is a breaking change for almost all functions, because the interface is now built around `CoupledSDEs`. The benefit is that the package now integrates much more seamlessly with DynamicalSystems.jl and DifferentialEquations.jl. To use the package with the old `StochSystem` struct, choose version v0.2.1 or lower.

Full changelog [here](https://github.com/JuliaDynamics/CriticalTransitions.jl/compare/v0.2.1...v0.3.0)

## v0.2.1
Freeze before major overhaul `v0.3.0`

#### Fixes
* bug fixes in the following functions:
    * `min_action_method`
    * `geometric_min_action_method`
    * `StochSystem(::CouplesODEs)`
    * `fhn_pathspace_sampling`

#### Enhancements
* expanded output of `transitions` function

#### Additions
* tests for various functions
* dev functions `basboundary`, `bisect_to_edge2`, `langevinmcmc_not_every_step`

#### Deprecations
* `to_cds` function (still exists but will raise deprecation warning)

Since `v0.2.0`, we also joined
[JuliaDynamics](https://juliadynamics.github.io/JuliaDynamics/) with this package, such that
the code is now hosted at [https://github.com/JuliaDynamics/CriticalTransitions.jl](https://github.com/JuliaDynamics/CriticalTransitions.jl), and improved the documentation.

## v0.2.0
Compatibility upgrade to `DynamicalSystems v3`

#### BREAKING changes
* changed `StochSystem.dim` to `StochSystem.u`
* changed keyword arguments in `basins` and `bisect_to_edge` functions
* removed `tocds` function (replaced by `CoupledODEs`)
* renamed functions:
    * `mam` -> `min_action_method`
    * `gmam` -> `geometric_min_action_method`
    * `FitzHughNagumo` -> `fitzhugh_nagumo`
    * `FitzHughNagumo!` -> `fitzhugh_nagumo!`

#### New functions
* `stochastic_bridge`
* `CoupledODEs(sys::StochSystem, ...)`
* `StochSystem(ode::CoupledODEs, ...)`

First prototype of a `RateSystem`.

## v0.1.1
Feature freeze before v0.2.0

This version is compatible with `v2` of `DynamicalSystems.jl` (specifically `DynamicalSystems v2.3.2`). The next release (`v0.2.0`) will include some breaking changes in the system type definitions and upgrade its dependencies to `DynamicalSystems v3`.

#### New functionality since previous release
* Large deviation theory: `fw_action`, `om_action`, `geometric_action`, `action`, `mam` and `gmam`
* Langevin MCMC sampling in pathspace: `langevinmcmc`
* Edge tracking algorithm: `edgetracking`
* Basins of attraction: `basins`, `basinboundary`
* additional convenience functions added

#### New predefined systems
* Ocean models: `stommel`, `cessi`, `rooth_smooth`
* Population dynamics: `originaltruscottbrindley`, `modifiedtruscottbrindley`, `rivals`

> Tested main functions on the FitzHugh-Nagumo model by running `test/functest.jl`.

## v0.1.0
First release! 🎉
> Tested for out-of-place system function `FitzHughNagumo` with WhiteGauss noise

#### Current functionality
Introduced the `StochSystem` type with the following methods:

* `equilib` - get equilibrium starting from initial condition
* `fixedpoints` - get fixedpoints of system as in _DynamicalSystems.jl_
* `relax` - simulate deterministic dynamics from initial condition
* `simulate` - simulate stochastic dynamics from initial condition
* `transition` - simulate transition sample from state 1 to state 2
* `transitions` - simulate an ensemble of transition samples from state 1 to 2
* `tocds` - convert `StochSystem` to `ContinuousDynamicalSystem` of _DynamicalSystems.jl_

Added convenience functions `make_jld2`, `sys_string`, `sys_info`, `idfunc`, `idfunc!`.

#### Systems
* `FitzHughNagumo`, `FitzHughNagumo!`