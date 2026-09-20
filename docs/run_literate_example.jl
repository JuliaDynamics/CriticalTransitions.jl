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

function preprocess(content)
    return replace(
        content,
        r"^#note # (.*)$"m => s"""
            # !!! note
            #     \1""",
    )
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
