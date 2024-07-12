using Documenter
using DocumenterCitations
using DocumenterInterLinks
using Pkg

using CriticalTransitions, ChaosTools, Attractors
using CairoMakie

project_toml = Pkg.TOML.parsefile(joinpath(@__DIR__, "..", "Project.toml"))
package_version = project_toml["version"]
name = project_toml["name"]
authors = join(project_toml["authors"], ", ") * " and contributors"
github = "https://github.com/juliadynamics/CriticalTransitions.jl"

links = InterLinks(
    "DiffEqNoiseProcess" => "https://docs.sciml.ai/DiffEqNoiseProcess/stable/",
    "DifferentialEquations" => "https://docs.sciml.ai/DiffEqDocs/stable/",
    "StochasticDiffEq" => "https://docs.sciml.ai/DiffEqDocs/stable/",
)

bib = CitationBibliography(joinpath(@__DIR__, "src", "refs.bib"); style=:numeric)

include("pages.jl")

html_options = Dict(
    :prettyurls => true,
    :canonical => "https://juliadynamics.github.io/CriticalTransitions.jl/",
    :mathengine => Documenter.MathJax2(),
)

if Documenter.DOCUMENTER_VERSION >= v"1.3.0"
    html_options[:inventory_version] = package_version
end

makedocs(;
    authors=authors,
    sitename="CriticalTransitions.jl",
    linkcheck=true,
    modules=[
        CriticalTransitions,
        Base.get_extension(CriticalTransitions, :ChaosToolsExt),
        Base.get_extension(CriticalTransitions, :CoupledSDEsBaisin),
    ],
    doctest=false,
    format=Documenter.HTML(; html_options...),
    warnonly=[:doctest, :missing_docs, :cross_references, :linkcheck],
    pages=pages,
    plugins=[bib, links],
)

deploydocs(; repo="github.com/JuliaDynamics/CriticalTransitions.jl.git", push_preview=false)
