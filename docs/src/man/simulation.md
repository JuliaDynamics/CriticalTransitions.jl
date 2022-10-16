# Simulating the system

## Features
Currently the following functions are implemented:
* `equilib(sys, state; kwargs...)`: returns the stable equilibrium of the initial point `state`.
* `fixedpoints(sys, box)`: wrapper for `fixedpoints` function of DynamicalSystems
* `simulate(sys, init; kwargs...)`: integrates the stochastic system forward in time
* `relax(sys, init; kwargs...)`: integrates the deterministic part of the system forward in time (evolution in the absence of noise)
* `transition(sys, x_i, x_f; kwargs...)`: simulate a sample transition trajectory from point x_i to point x_f
* `transitions(sys, x_i, x_f, N=1; kwargs)`: simulate an ensemble of N sample transitions
