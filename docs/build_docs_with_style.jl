CI = get(ENV, "CI", nothing) == "true" || get(ENV, "GITHUB_TOKEN", nothing) !== nothing

import Pkg
Pkg.pkg"add Documenter@1"

# Load documenter
using Documenter
using DocumenterTools: Themes
ENV["JULIA_DEBUG"] = "Documenter"

# For easier debugging when downloading from a specific branch.
github_user = "JuliaDynamics"
branch = "master"
download_path = "https://raw.githubusercontent.com/$github_user/doctheme/$branch"

import Downloads
for file in ("juliadynamics-lightdefs.scss", "juliadynamics-darkdefs.scss", "juliadynamics-style.scss")
    Downloads.download("$download_path/$file", joinpath(@__DIR__, file))
end

# create the themes
for w in ("light", "dark")
    header = read(joinpath(@__DIR__, "juliadynamics-style.scss"), String)
    theme = read(joinpath(@__DIR__, "juliadynamics-$(w)defs.scss"), String)
    write(joinpath(@__DIR__, "juliadynamics-$(w).scss"), header * "\n" * theme)
end

# compile the themes
Themes.compile(joinpath(@__DIR__, "juliadynamics-light.scss"), joinpath(@__DIR__, "src/assets/themes/documenter-light.css"))
Themes.compile(joinpath(@__DIR__, "juliadynamics-dark.scss"), joinpath(@__DIR__, "src/assets/themes/documenter-dark.css"))

# Download and apply CairoMakie plotting style
using CairoMakie
Downloads.download("$download_path/style.jl", joinpath(@__DIR__, "style.jl"))
include("style.jl")

"""
    build_docs_with_style(pages::Vector, modules... ; kw...)

Call the `makedocs` function with some predefined style components.
The first module dictates site name, while the rest need to be included
to expand and cross-referrence docstrings from other modules.
`kw` are propagated to `makedocs` while the keyword `htmlkw` is propagated to
`Documenter.HTML`.
"""
function build_docs_with_style(
        pages, modules...;
        bib = nothing, plugins = nothing, authors = "George Datseris and contributors",
        htmlkw = NamedTuple(), kwargs...
    )

    if isnothing(plugins)
        if !isnothing(bib)
            plugins = [bib,]
        else
            plugins = []
        end
    end

    makedocs(;
        modules = [modules...],
        format = Documenter.HTML(;
            prettyurls = CI,
            assets = [
                asset("https://fonts.googleapis.com/css?family=Montserrat|Source+Code+Pro&display=swap", class = :css),
            ],
            collapselevel = 3,
            htmlkw...,
        ),
        sitename = "$(modules[1]).jl",
        authors,
        pages,
        checkdocs = :exported,
        linkcheck = true,
        linkcheck_timeout = 2,
        expandfirst = ["index.md"],
        doctest = false,
        # The following Documenter fails will NOT ERROR the docbuild!
        warnonly = [:doctest],
        kwargs..., plugins,
    )

    return if CI
        deploydocs(
            repo = "github.com/JuliaDynamics/$(modules[1]).jl.git",
            target = "build",
            push_preview = true
        )
    end
end
