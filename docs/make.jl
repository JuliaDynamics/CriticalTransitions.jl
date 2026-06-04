# Loading
CI = get(ENV, "CI", nothing) == "true" || get(ENV, "GITHUB_TOKEN", nothing) !== nothing
using CriticalTransitions, ChaosTools, Attractors
using DynamicalSystemsBase
using StochasticDiffEq, DiffEqNoiseProcess

# Literate-render `examples/*.jl` into `src/examples/*.md` on every build so
# local previews stay in sync with CI. The generated `.md` files are gitignored.
include("make_md_examples.jl")

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
    "ChaosTools" => "https://juliadynamics.github.io/ChaosTools.jl/stable/",
    # "DiffEqNoiseProcess" => "https://docs.sciml.ai/DiffEqNoiseProcess/stable/",
    # "DifferentialEquations" => "https://docs.sciml.ai/DiffEqDocs/stable/",
    # "StochasticDiffEq" => "https://docs.sciml.ai/DiffEqDocs/stable/",
)

# `remotes` tells Documenter where the source of each included module lives on
# GitHub, so the "source" link on every docstring resolves to the correct
# repository and version tag. We include DynamicalSystemsBase and Attractors in
# `modules` below (their docstrings are rendered here), so without `remotes`
# Documenter would either fail or link to a wrong path. The version string is
# read from `Pkg.status` so the link points at the actually installed version.
using Pkg
project_toml = Pkg.TOML.parsefile(joinpath(@__DIR__, "..", "Project.toml"))
package_version = project_toml["version"]
authors = join(project_toml["authors"], ", ") * " and contributors"

remote_pairs = map(["DynamicalSystemsBase", "Attractors", "ChaosTools"]) do pkg_name
    status = sprint(io -> Pkg.status(pkg_name; io = io))
    version = match(r"(v[0-9].[0-9]+.[0-9]+)", status)[1]
    gh_moi = Documenter.Remotes.GitHub("JuliaDynamics", string(pkg_name, ".jl"))
    (pkgdir(eval(Symbol(pkg_name))) => (gh_moi, version))
end
remotes = Dict(remote_pairs)

# Extra kwargs forwarded to `Documenter.HTML(...)` inside `build_docs_with_style`.
# The doctheme already sets `prettyurls = CI`, `collapselevel = 3`, and assets,
# so only add things that genuinely differ from those defaults:
#   canonical  -> SEO canonical URL for the deployed `stable/` site.
#   mathengine -> MathJax 2 instead of Documenter's default KaTeX. The current
#                 math in `src/` (\mathbb, \mathcal, aligned, ...) is supported
#                 by both engines; MathJax is kept as a safety net for richer
#                 LaTeX (\bm, \tag, AMS environments) we may add later.
html_options = Dict(
    :canonical => "https://juliadynamics.github.io/CriticalTransitions.jl/",
    :mathengine => Documenter.MathJax2(),
    # :example_size_threshold => nothing,
    # :size_threshold_warn => nothing,
    # :size_threshold => nothing,
)

# `inventory_version` tags the `objects.inv` file Documenter generates so that
# downstream docs sites using DocumenterInterLinks can pin cross-references to
# a specific CriticalTransitions release. The field only exists in Documenter
# >= 1.3, hence the version guard.
if Documenter.DOCUMENTER_VERSION >= v"1.3.0"
    html_options[:inventory_version] = package_version
end

modules = [CriticalTransitions, DynamicalSystemsBase, Attractors, ChaosTools]

# URLs that are valid but cannot pass `linkcheck`: academic publishers that
# return 403 to automated requests (Wiley, Annual Reviews) and a few external
# pages that are reachable in a browser but unreliable under the short
# `linkcheck_timeout` (author homepages, the CriticalEarth site, the SciML
# DiffEq redirect). Skipping the check for these keeps the build deterministic.
linkcheck_ignore = [
    r"https://onlinelibrary\.wiley\.com/.*",
    r"https://www\.annualreviews\.org/.*",
    r"https://homepages\.warwick\.ac\.uk/.*",
    r"https://www\.math\.drexel\.edu/.*",
    r"https://diffeq\.sciml\.ai/.*",
    r"https://www\.criticalearth\.eu.*",
]

# Build docs
build_docs_with_style(
    pages,
    modules...;
    plugins = [bib, links],
    authors,
    htmlkw = html_options,
    checkdocs_ignored_modules = [CriticalTransitions.CTLibrary],
    linkcheck_ignore,
    remotes,
)
