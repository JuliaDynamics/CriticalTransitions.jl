using BenchmarkTools
using CriticalTransitions

const SUITE = BenchmarkGroup()

include("kpo.jl")
include("maierstein.jl")
include("multiplicative_noise.jl")
include("gmam_lqa.jl")
# Disabled until #280's RateSystem refactor is finished — see PR #310.
# include("ratesystem.jl")

benchmark_KPO!(SUITE)
benchmark_maierstein!(SUITE)
benchmark_multiplicative_noise!(SUITE)
benchmark_gmam_lqa!(SUITE)
# benchmark_rate_system!(SUITE)

BenchmarkTools.tune!(SUITE)
results = BenchmarkTools.run(SUITE; verbose = true)
display(median(results))

BenchmarkTools.save("benchmarks_output.json", median(results))
