using Documenter, CriticalTransitions

makedocs(sitename="CriticalTransitions.jl", doctest=false)

deploydocs(
    repo = "github.com/reykboerner/CriticalTransitions.jl.git",
)
