using Literate

length(ARGS) == 2 || error("usage: julia run_literate_example.jl INPUT.jl OUTPUT_DIR")

input_file = abspath(ARGS[1])
output_dir = abspath(ARGS[2])
repo_root = normpath(joinpath(@__DIR__, ".."))

extra_literate_config = if isempty(get(ENV, "CI", ""))
    Dict("repo_root_path" => repo_root, "repo_root_url" => "file://" * repo_root)
else
    Dict()
end

const PNG_WRAPPER = """
struct _LiteratePNGOnly{T} #hide
    value::T #hide
end #hide
Base.showable(::MIME\"text/html\", ::_LiteratePNGOnly) = false #hide
Base.showable(::MIME\"image/svg+xml\", ::_LiteratePNGOnly) = false #hide
Base.showable(::MIME\"text/markdown\", ::_LiteratePNGOnly) = false #hide
Base.showable(::MIME\"image/png\", x::_LiteratePNGOnly) = showable(MIME(\"image/png\"), x.value) #hide
Base.show(io::IO, mime::MIME\"image/png\", x::_LiteratePNGOnly) = show(io, mime, x.value) #hide
"""

function preprocess(content)
    # Preserve the package's existing `#note` shorthand.
    content = replace(
        content,
        r"^#note # (.*)$"m => s"""
            # !!! note
            #     \1""",
    )

    # With execute=true, DocumenterFlavor would otherwise prefer Makie's
    # text/html representation. Keep the visible `fig` expression unchanged,
    # but make the hidden final expression expose only image/png so generated
    # pages stay compact and retain Documenter @ref/@id semantics.
    content = replace(
        content,
        r"(?m)^(fig[A-Za-z0-9_]*)[ \t]*$" => s"\1\n_LiteratePNGOnly(\1) #hide",
    )
    return PNG_WRAPPER * "\n" * content
end

Literate.markdown(
    input_file,
    output_dir;
    flavor = Literate.DocumenterFlavor(),
    credit = false,
    config = extra_literate_config,
    execute = true,
    preprocess,
)
