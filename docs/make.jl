# Loading
CI = get(ENV, "CI", nothing) == "true" || get(ENV, "GITHUB_TOKEN", nothing) !== nothing
using CriticalTransitions, ChaosTools, Attractors
using DynamicalSystemsBase
using StochasticDiffEq, DiffEqNoiseProcess

# TODO: It is bad practice to access extensions at the top level.
# We need to remove these, but not clear where they are used currently.
StochasticSystemsBase = Base.get_extension(DynamicalSystemsBase, :StochasticSystemsBase)
ChaosToolsExt = Base.get_extension(CriticalTransitions, :ChaosToolsExt)
CoupledSDEsBasin = Base.get_extension(CriticalTransitions, :CoupledSDEsBasin)

if get(ENV, "CI", nothing) == "true"
    include("make_md_examples.jl")
end

include("pages.jl")

# Install style  of JuliaDynamics
using Downloads: Downloads
Downloads.download(
    "https://raw.githubusercontent.com/JuliaDynamics/doctheme/master/build_docs_with_style.jl",
    joinpath(@__DIR__, "build_docs_with_style.jl"),
)
include("build_docs_with_style.jl")

# Set up doc build options
using DocumenterCitations
using DocumenterInterLinks

bib = CitationBibliography(joinpath(@__DIR__, "src", "refs.bib"); style = :authoryear)

links = InterLinks(
    "DynamicalSystemsBase" => "https://juliadynamics.github.io/DynamicalSystemsBase.jl/stable/",
    "Attractors" => "https://juliadynamics.github.io/Attractors.jl/stable/",
    # "DiffEqNoiseProcess" => "https://docs.sciml.ai/DiffEqNoiseProcess/stable/",
    # "DifferentialEquations" => "https://docs.sciml.ai/DiffEqDocs/stable/",
    # "StochasticDiffEq" => "https://docs.sciml.ai/DiffEqDocs/stable/",
)

# TODO: What is this? What is this remotes? Why do we need it? Seems like really complex code.
using Pkg
project_toml = Pkg.TOML.parsefile(joinpath(@__DIR__, "..", "Project.toml"))
package_version = project_toml["version"]
authors = join(project_toml["authors"], ", ") * " and contributors"

remote_pairs = map(["DynamicalSystemsBase", "Attractors"]) do pkg_name
    status = sprint(io -> Pkg.status(pkg_name; io = io))
    version = match(r"(v[0-9].[0-9]+.[0-9]+)", status)[1]
    gh_moi = Documenter.Remotes.GitHub("JuliaDynamics", string(pkg_name, ".jl"))
    (pkgdir(eval(Symbol(pkg_name))) => (gh_moi, version))
end
remotes = Dict(remote_pairs)

# TODO: Why do we need this?
html_options = Dict(
    :prettyurls => true,
    :canonical => "https://juliadynamics.github.io/CriticalTransitions.jl/",
    :mathengine => Documenter.MathJax2(),
    # :example_size_threshold => nothing,
    # :size_threshold_warn => nothing,
    # :size_threshold => nothing,
)

# TODO: Why do we need these?
if Documenter.DOCUMENTER_VERSION >= v"1.3.0"
    html_options[:inventory_version] = package_version
end

modules = filter(
    !isnothing,
    [
        CriticalTransitions,
        DynamicalSystemsBase,
        Attractors,
        Base.get_extension(CriticalTransitions, :ChaosToolsExt),
        Base.get_extension(CriticalTransitions, :AttractorsExt),
        Base.get_extension(DynamicalSystemsBase, :StochasticSystemsBase),
    ],
)

# Build docs
build_docs_with_style(
    pages,
    modules...;
    plugins = [bib, links],
    # Broad warnonly until #280's docs cleanup lands; some example pages have
    # stale code that errors at @example evaluation.
    warnonly = true,
    authors,
    htmlkw = html_options,
    checkdocs_ignored_modules = [CriticalTransitions.CTLibrary],
    remotes,
)
