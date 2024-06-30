push!(LOAD_PATH, "../src/")

using Documenter
using DocumenterCitations

using CriticalTransitions, ChaosTools, Attractors
using CairoMakie

bib = CitationBibliography(joinpath(@__DIR__, "src", "refs.bib"); style=:numeric)
include("pages.jl")

makedocs(;
    sitename="CriticalTransitions.jl",
    repo=Documenter.Remotes.GitHub("JuliaDynamics", "CriticalTransitions.jl"),
    modules=[
        CriticalTransitions,
        Base.get_extension(CriticalTransitions, :ChaosToolsExt),
        Base.get_extension(CriticalTransitions, :CoupledSDEsBaisin),
    ],
    doctest=false,
    format=Documenter.HTML(;
        canonical  = "https://juliadynamics.github.io/CriticalTransitions.jl/",
        prettyurls = true,
        mathengine = Documenter.MathJax2(),
    ),
    linkcheck=true,
    warnonly=[:doctest, :missing_docs, :cross_references, :linkcheck],
    pages=pages,
    plugins=[bib],
)

deploydocs(; repo="github.com/JuliaDynamics/CriticalTransitions.jl.git", push_preview=false)
