# CriticalTransitions.jl

```@docs
CriticalTransitions
```

## Installation
!!! note "Installing Julia"
    To install the Julia language, we recommend [juliaup](https://github.com/JuliaLang/juliaup).

`CriticalTransitions` is a registered Julia package and can be installed with the Julia package manager:

```console
julia> ]
Pkg>   add CriticalTransitions
```
or 

```julia
using Pkg; Pkg.add("CriticalTransitions")
```

## Getting started

The [Tutorial](@ref "CriticalTransitions.jl Tutorial") showcases how to set up different forced dynamical systems and how to analyze their transition behavior using some of the key functionality. For more detailed and advanced examples, see the **Examples** section. The **Manual** section lists the full API with additional explanation.

## List of methods

```@index
Pages = map(file -> joinpath("man", file), readdir("man"))
Order   = [:type, :function]
```

## People

Main developers:
- Reyk Börner (@reykboerner)
- Orjan Ameye (@oameye)
- Ryan Deeley (@ryandeeley)
- Raphael Römer (@raphael-roemer)
- George Datseris (@Datseris)

Thanks to Jeroen Wouters, Calvin Nesbitt, Tobias Grafke & Oliver Mehling.

This package got started in the EU-funded [CriticalEarth](https://www.criticalearth.eu) project.