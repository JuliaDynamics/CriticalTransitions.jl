using Literate

### Process examples
# Always rerun examples
const EXAMPLES_IN = joinpath(@__DIR__, "..", "examples")
const OUTPUT_MD_DIR = joinpath(@__DIR__, "src", "examples")
const LITERATE_RUNNER = joinpath(@__DIR__, "run_literate_example.jl")
const DOCS_PROJECT = abspath(@__DIR__)

# Keep plotting backends headless in both CI and local documentation builds.
get!(ENV, "GKSwstype", "100")

examples = filter!(file -> endswith(file, ".jl"), readdir(EXAMPLES_IN; join = true))
filter!(file -> !contains(file, "make_nb_examples"), examples)
# Orphan / WIP sources kept in `examples/` for reference but not part of the
# rendered docs. Skip them so Documenter doesn't process broken Literate output.
const SKIP_EXAMPLES = ("Langevin_MCMC",)
filter!(file -> !any(name -> occursin(name, file), SKIP_EXAMPLES), examples)
sort!(examples)
const EXAMPLE_STEMS = Set(first(splitext(basename(example))) for example in examples)

const LITERATE_WORKERS = let
    workers = try
        parse(Int, get(ENV, "DOCUMENTER_LITERATE_WORKERS", "4"))
    catch
        error("DOCUMENTER_LITERATE_WORKERS must be a positive integer")
    end
    workers > 0 || error("DOCUMENTER_LITERATE_WORKERS must be a positive integer")
    workers
end

function generated_output(name)
    return any(
        startswith(name, stem * ".") || startswith(name, stem * "-") for
            stem in EXAMPLE_STEMS
    )
end

# Remove only artifacts belonging to the pages generated below. This also
# clears images left behind when an example changes its output blocks.
mkpath(OUTPUT_MD_DIR)
for file in readdir(OUTPUT_MD_DIR)
    generated_output(file) && rm(joinpath(OUTPUT_MD_DIR, file); force = true)
end

# Each task owns a separate Julia process. This keeps example state isolated and
# gives bounded process-level parallelism without sharing plotting/output state.
function generate_example(example)
    run(
        `$(Base.julia_cmd()) --startup-file=no --project=$(DOCS_PROJECT) --threads=1 $(LITERATE_RUNNER) $(example) $(OUTPUT_MD_DIR)`,
    )
    return example
end

if !isempty(examples)
    asyncmap(generate_example, examples; ntasks = min(LITERATE_WORKERS, length(examples)))
end
