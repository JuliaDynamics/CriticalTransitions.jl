# CriticalTransitions.jl

A Julia package for the numerical investigation of critical transitions in metastable systems, building on `DifferentialEquations.jl` and `DynamicalSystems.jl`.

## Usage
As this module is not published yet, it needs to be accessed locally (rather than using `Pkg`.)

### Getting started
Clone the repo using `git clone https://github.com/reykboerner/CriticalTransitions.jl.git`.

### Loading the module
1. In your Julia script, include the CriticalTransitions source code: `include("PATH/src/CriticalTransitions.jl")`. Here `PATH` is your local relative path to the `CriticalTransitions.jl` folder.
2. Load the module: `using .CriticalTransitions`.

### Setting up a StochSystem
An instance of a stochastic dynamical system is defined via `StochSystem(f, pf, dim [, σ [, g, pg, Σ , process]])`,
taking the following arguments:
* `f` (Function): Dynamical ODE rule describing the drift term of the system, corresponding to `f` in the ODEProblem of `DifferentialEquations`. Can be defined in-place (`f!(du, u, p, t)`) or out-of-place (`f(u,p,t)`).
* `pf` (Vector or Nothing): Parameter vector for the drift term.
* `dim` (Int64): Dimension of the system.
* `σ` (Float64): Noise intensity. Defaults to `1.0`.
* `g` (Function): Dynamical ODE rule describing the noise term of the system. Same format as `f`. Defaults to `idfunc()`.
* `pg` (Vector or Nothing): Parameter vector of the noise term.
* `Σ` (Matrix): Noise covariance matrix. Defaults to `I` (identity matrix).
* `process` (String): Noise process. Defaults to `white-gauss` (independent n-dimensional Brownian motion).

> Note: Some dynamical systems are predefined in the `systems` folder and can be inserted for `f`.

### Analyzing StochSystems
Currently the following functions are implemented:
* `equilib(sys, state; kwargs...)`: returns the stable equilibrium of the initial point `state`.
* `fixedpoints(sys, box)`: wrapper for `fixedpoints` function of DynamicalSystems
* `simulate(sys, init; kwargs...)`: integrates the stochastic system forward in time
* `relax(sys, init; kwargs...)`: integrates the deterministic part of the system forward in time (evolution in the absence of noise)

### Example: Bistable FitzHugh-Nagumo model
```
include("src/CriticalTransitions.jl")
include("systems/fitzhughnagumo.jl")

using .CriticalTransitions

# Parameters
p = [1., 3., 1., 1., 1., 0.]

# StochSystem
sys = StochSystem(fitzhughnagumo, p, 2)

# Find equilibrium starting from point [1,1]
A = equilib(sys, [1,1])
```

### Git cheat sheet

remote = the GitHub server; local = your computer

#### Before any editing session
* `git remote update` to refresh info about changes on remote
* `git status` to check whether you are ahead of, behind, or in sync with remote
* `git pull` to get the latest version from remote

> Never use `git pull` if there are still uncommitted local changes!

#### Checking what you've changed locally
* `git status` to see status of locally modified files
* `git diff` to see what the changes actually are line by line within the files

> If you just want to see the changes within a specific file, type `git diff <filename>`

#### Ready to add, commit & push
* `git add .` to add all unstaged changes (or `git add <file1> <file2> ...` to only add certain files)
* `git commit -m "<message>"` to commit the added changes
* `git push` to push to remote (github)