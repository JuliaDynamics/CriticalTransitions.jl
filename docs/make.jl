using Documenter
using DocumenterCitations
using CriticalTransitions

bib = CitationBibliography(
    joinpath(@__DIR__, "src", "refs.bib");
    style=:numeric
)

makedocs(;
    sitename="CriticalTransitions.jl",
    repo = Documenter.Remotes.GitHub("JuliaDynamics", "CriticalTransitions.jl"),
    modules=[CriticalTransitions],
    doctest=false,
    format = Documenter.HTML(
        canonical="https://juliadynamics.github.io/CriticalTransitions.jl/",
        prettyurls = true,
        mathengine = Documenter.MathJax2()
        ),
    linkcheck = true,
    warnonly = [:doctest, :missing_docs, :cross_references, :linkcheck],
    pages=Any[
        "Home" => "index.md",
        "Quickstart" => "quickstart.md",
        "Tutorial" => "examples/tutorial.md",
        "Manual" => Any[
            "Defining a StochSystem" => "man/stochsystem.md",
            "Stability analysis" => "man/systemanalysis.md",
            "Simulating the system" => "man/simulation.md",
            "Sampling transitions" => "man/sampling.md",
            "Large deviation theory" => "man/largedeviations.md",
            "Noise processes" => "man/noise.md",
            "Utilities" => "man/utils.md"
        ],
        "Examples" => Any[
            "Instanton for the Maier-Stein system" => "examples/gMAM_Maierstein.md"
        ],
        "Predefined systems" => "man/systems.md",
        "Development stage" => "man/dev.md"
    ],
    plugins=[bib]
)

deploydocs(
    repo = "github.com/JuliaDynamics/CriticalTransitions.jl.git",
    push_preview = false
)
