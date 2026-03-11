# CriticalTransitions.jl

[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://juliadynamics.github.io/CriticalTransitions.jl/dev/)
[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://juliadynamics.github.io/CriticalTransitions.jl/stable/)
[![Tests](https://github.com/JuliaDynamics/CriticalTransitions.jl/actions/workflows/Tests.yml/badge.svg)](github.com/JuliaDynamics/CriticalTransitions.jl/actions/workflows/ci.yml)
[![Benchmarks](https://github.com/JuliaDynamics/CriticalTransitions.jl/actions/workflows/Benchmarks.yaml/badge.svg?branch=main)](https://juliadynamics.github.io/CriticalTransitions.jl/benchmarks/)

[![Code Style: Blue](https://img.shields.io/badge/code%20style-blue-4495d1.svg)](https://github.com/invenia/BlueStyle)
[![Aqua QA](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)
[![JET](https://img.shields.io/badge/%E2%9C%88%EF%B8%8F%20tested%20with%20-%20JET.jl%20-%20red)](https://github.com/aviatesk/JET.jl)

**CriticalTransitions.jl** is a first-of-its-kind software for formalizing, automating, and making extendable, the analysis of critical transitions in dynamical systems.
Current content exists along two independent paths: noise- and rate- induced transitions.

The main software highlights are:

- easily construct stochastic and nonautonomous dynamical systems
- efficiently sample transition path ensembles
- calculate minimum action paths and critical forcing rates
- use a growing toolbox of tested and documented functions implementing concepts
  of large deviation theory, transition path theory, and rate-induced tipping
- and more features shown in the documentation and planned for the future!

CriticalTransitions.jl can be used as a standalone package, or as part of
[DynamicalSystems.jl](https://juliadynamics.github.io/DynamicalSystemsDocs.jl/dynamicalsystems/stable/).

To install it, run `import Pkg; Pkg.add("CriticalTransitions")`.

All further information is provided in the documentation, which you can either find [online](https://juliadynamics.github.io/DynamicalSystemsDocs.jl/criticaltransitions/stable/) or build locally by running the `docs/make.jl` file.
