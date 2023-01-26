using Documenter, CriticalTransitions

makedocs(
    sitename="CriticalTransitions.jl",
    doctest=false,
    format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true",
        mathengine = Documenter.MathJax2()),
    
    pages=Any[
        "Home" => "index.md",
        "Quickstart" => "quickstart.md",
        "Tutorial" => "tutorial.md",
        "Manual" => Any[
            "Defining a StochSystem" => "man/stochsystem.md",
            "Stability analysis" => "man/systemanalysis.md",
            "Simulating the system" => "man/simulation.md",
            "Sampling transitions" => "man/sampling.md",
            "Large deviation theory" => "man/largedeviations.md",
            "Noise processes" => "man/noise.md",
            "Utilities" => "man/utils.md"
        ],
        "Predefined systems" => "man/systems.md",
        "Development stage" => "man/dev.md"
    ]
)

deploydocs(
    repo = "github.com/reykboerner/CriticalTransitions.jl.git",
)