# CriticalTransitions.jl

A Julia package for the numerical investigation of critical transitions in metastable systems, builiding upon `DifferentialEquations.jl` and `DynamicalSystems.jl`.

## Usage
As this module is not published yet, it needs to be accessed locally (rather than using `Pkg`.)

### Getting started
Clone the repo using `git clone https://github.com/reykboerner/CriticalTransitions.jl.git`.

### Loading the module
1. In your Julia script, include the CriticalTransitions source code: `include("PATH/src/CriticalTransitions.jl")`. Here `PATH` is your local path to the `CriticalTransitions.jl` folder.
2. Load the module: `using .CriticalTransitions`.

### Setting up a StochSystem
An instance of a stochastic dynamical system is defined via `StochSystem(f, p, [$\sigma$, g, pg, $\Sigma$, process])`,
taking the following arguments:
* `f` (Function): Dynamical ODE rule describing the drift term of the system, corresponding to `f` in the ODEProblem of `DifferentialEquations`. Can be defined in-place (`f!(du, u, p, t)`) or out-of-place (`f(u,p,t)`).
* `p` (Vector or Nothing): Parameter vector for the drift term.
* `$\sigma$` (Float64): Noise intensity. Defaults to `1.0`.
* `g` (Function): Dynamical ODE rule describing the noise term of the system. Same format as `f`. Defaults to `idfunc()`.
* `pg` (Vector or Nothing): Parameter vector of the noise term.
* `$\Sigma$` (Matrix): Noise covariance matrix. Defaults to `I` (identity matrix).
* `process` (String): Noise process. Defaults to `white-gauss` (independent n-dimensional Brownian motion).

> Note: Some dynamical systems are predefined in the `systems` folder and can be inserted for `f`.

### Analyzing StochSystems
Currently the following functions are implemented:
* `equilib(sys, state; kwargs...)`: returns the stable equilibrium of the initial point `state`.

### Example: Bistable FitzHugh-Nagumo model
```
include("src/CriticalTransitions.jl")
include("systems/fitzhughnagumo.jl")

using .CriticalTransitions

# Parameters
p = [1., 3., 1., 1., 1., 0.]

# StochSystem
sys = StochSystem(fitzhughnagumo, p)

# Find equilibrium starting from point [1,1]
A = equilib(sys, [1,1])
```
