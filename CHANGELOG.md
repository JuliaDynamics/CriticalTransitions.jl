# Changelog for `CriticalTransitions.jl`

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