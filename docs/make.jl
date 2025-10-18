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

remote_pairs = map(["DynamicalSystemsBase", "Attractors"]) do pkg_name
    status = sprint(io -> Pkg.status(pkg_name; io=io))
    version = match(r"(v[0-9].[0-9]+.[0-9]+)", status)[1]
    gh_moi = Documenter.Remotes.GitHub("JuliaDynamics", string(pkg_name, ".jl"))
    (pkgdir(eval(Symbol(pkg_name))) => (gh_moi, version))
end
remotes = Dict(remote_pairs)

makedocs(;
    authors,
    sitename="CriticalTransitions.jl",
    repo=Documenter.Remotes.GitHub("JuliaDynamics", "CriticalTransitions.jl"),
    modules=[
        CriticalTransitions,
        DynamicalSystemsBase,
        Attractors,
        Base.get_extension(CriticalTransitions, :ChaosToolsExt),
        Base.get_extension(CriticalTransitions, :AttractorsExt),
        Base.get_extension(DynamicalSystemsBase, :StochasticSystemsBase),
    ],
    remotes,
    checkdocs_ignored_modules=[CriticalTransitions.CTLibrary],
    pages,
    linkcheck=true,
    pagesonly=true,
    checkdocs=:exports,
    doctest=false,
    format=Documenter.HTML(; html_options...),
    warnonly=[:missing_docs, :linkcheck, :cross_references],
    plugins=[bib, links],
    draft=CI ? false : true,
)

deploydocs(; repo="github.com/JuliaDynamics/CriticalTransitions.jl.git", push_preview=true)
