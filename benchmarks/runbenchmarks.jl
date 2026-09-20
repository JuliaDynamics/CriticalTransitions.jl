using BenchmarkTools
using CriticalTransitions

const SUITE = BenchmarkGroup()

include("transition_callbacks.jl")
benchmark_transition_callbacks!(SUITE)

BenchmarkTools.tune!(SUITE)
results = BenchmarkTools.run(SUITE; verbose = true)
display(median(results))

BenchmarkTools.save("benchmarks_output.json", median(results))
