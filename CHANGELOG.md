# Changelog for `CriticalTransitions.jl`

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