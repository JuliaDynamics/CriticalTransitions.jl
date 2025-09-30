CI = get(ENV, "CI", nothing) == "true" || get(ENV, "GITHUB_TOKEN", nothing) !== nothing

using Documenter
using DocumenterCitations
using DocumenterInterLinks
using Pkg

using CriticalTransitions, ChaosTools, Attractors
using DynamicalSystemsBase
StochasticSystemsBase = Base.get_extension(DynamicalSystemsBase, :StochasticSystemsBase)
ChaosToolsExt = Base.get_extension(CriticalTransitions, :ChaosToolsExt)
CoupledSDEsBasin = Base.get_extension(CriticalTransitions, :CoupledSDEsBasin)

using StochasticDiffEq, DiffEqNoiseProcess

project_toml = Pkg.TOML.parsefile(joinpath(@__DIR__, "..", "Project.toml"))
package_version = project_toml["version"]
name = project_toml["name"]
authors = join(project_toml["authors"], ", ") * " and contributors"
github = "https://github.com/juliadynamics/CriticalTransitions.jl"

links = InterLinks(
    "DynamicalSystemsBase" => "https://juliadynamics.github.io/DynamicalSystemsBase.jl/stable/",
    "Attractors" => "https://juliadynamics.github.io/Attractors.jl/stable/",
    # "DiffEqNoiseProcess" => "https://docs.sciml.ai/DiffEqNoiseProcess/stable/",
    # "DifferentialEquations" => "https://docs.sciml.ai/DiffEqDocs/stable/",
    # "StochasticDiffEq" => "https://docs.sciml.ai/DiffEqDocs/stable/",
)

bib = CitationBibliography(joinpath(@__DIR__, "src", "refs.bib"); style=:numeric)

if CI
    include("make_md_examples.jl")
else
    nothing
end

include("pages.jl")

html_options = Dict(
    :prettyurls => true,
    :canonical => "https://juliadynamics.github.io/CriticalTransitions.jl/",
    :mathengine => Documenter.MathJax2(),
    # :example_size_threshold => nothing,
    # :size_threshold_warn => nothing,
    # :size_threshold => nothing,
)

if Documenter.DOCUMENTER_VERSION >= v"1.3.0"
    html_options[:inventory_version] = package_version
end

makedocs(;
    debug=true,
    authors=authors,
    sitename="CriticalTransitions.jl",
    repo=Documenter.Remotes.GitHub("JuliaDynamics", "CriticalTransitions.jl"),
    modules=[
        CriticalTransitions,
        DynamicalSystemsBase,
        Attractors,
        Base.get_extension(CriticalTransitions, :ChaosToolsExt),
        Base.get_extension(CriticalTransitions, :AttractorsExt),
        Base.get_extension(DynamicalSystemsBase, :StochasticSystemsBase),
        # DynamicalSystemsBase
    ],
    checkdocs_ignored_modules=[CriticalTransitions.CTLibrary],
    pages=pages,
    linkcheck=true,
    #linkcheck_ignore = [r"http://docs\.juliadiffeq\.org/.*"],
    pagesonly=true,
    checkdocs=:exports,
    doctest=false,
    format=Documenter.HTML(; html_options...),
    warnonly=[:missing_docs, :linkcheck, :cross_references],
    plugins=[bib, links],
    # plugins=[links],
    draft=CI ? false : true,
)

deploydocs(; repo="github.com/JuliaDynamics/CriticalTransitions.jl.git", push_preview=true)
