using CriticalTransitions
using Test

using Random
const SEED = 0xd8e5d8df
Random.seed!(SEED)

files = [
    "gMAM.jl",
    "stochsystem.jl",
    "MAM.jl"
]

for file in files
    include(file)
    printstyled(file * ":    OK\n"; color = :green)
end
