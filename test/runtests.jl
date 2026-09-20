using CriticalTransitions
using ParallelTestRunner: ParallelTestRunner

# Start with autodiscovered tests and preserve the existing default suite.
testsuite = ParallelTestRunner.find_tests(@__DIR__)
args = ParallelTestRunner.parse_args(ARGS)

if ParallelTestRunner.filter_tests!(testsuite, args)
    # JET has its own dedicated workflow/environment, and API.jl is a legacy
    # file that was not part of the previous runtests.jl include list.
    delete!(testsuite, "quality/JET")
    delete!(testsuite, "largedeviations/API")
end

# Preserve the shared setup that the previous monolithic runtests.jl provided.
# ParallelTestRunner evaluates this in each isolated test sandbox.
const init_code = quote
    using CriticalTransitions, StaticArrays
    const CT = CriticalTransitions
    using Test
    using Random
    const SEED = 0xd8e5d8df
    Random.seed!(SEED)
    using CriticalTransitions.CTLibrary: fitzhugh_nagumo
end

ParallelTestRunner.runtests(CriticalTransitions, args; testsuite, init_code)
