# Simulating the system

We provide two main functions to simulate a [`CoupledSDEs`](@ref) forward in time:
* [`DynamicalSystemsBase.trajectory`](@extref), which integrates the stochastic `CoupledSDEs` system forward in time
* [`deterministic_orbit`](@ref), which integrates only the deterministic part of the `CoupledSDEs` system 


```@docs
deterministic_orbit
```