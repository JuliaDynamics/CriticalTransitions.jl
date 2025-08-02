using BenchmarkTools
using CriticalTransitions

const SUITE = BenchmarkGroup()

include("kpo.jl")
include("maierstein.jl")

benchmark_KPO!(SUITE)
benchmark_maierstein!(SUITE)

BenchmarkTools.tune!(SUITE)
results = BenchmarkTools.run(SUITE; verbose=true)
display(median(results))

BenchmarkTools.save("benchmarks_output.json", median(results))
