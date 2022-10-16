# CriticalTransitions.jl

[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://reykboerner.github.io/CriticalTransitions.jl/dev/)

A Julia package for the numerical investigation of critical transitions in metastable systems, building on `DifferentialEquations.jl` and `DynamicalSystems.jl`.

## Usage
See [package documentation](https://reykboerner.github.io/CriticalTransitions.jl/dev/).

## Example: Bistable FitzHugh-Nagumo model
```julia
include("src/CriticalTransitions.jl")
using .CriticalTransitions

# Parameters
p = [1., 3., 1., 1., 1., 0.]

# StochSystem
sys = StochSystem(FitzHughNagumo, p, 2)

# Get stable fixed points
eqs, eigs, stab = fixedpoints(sys, [-2,-2], [2,2])
fp1, fp2 = eqs[stab]

# Simulate noisy trajectory starting from fixed point 1
sim = simulate(sys, fp1)
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
