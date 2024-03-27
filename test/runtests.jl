using CriticalTransitions, StaticArrays
CT = CriticalTransitions
using Test

using Random
const SEED = 0xd8e5d8df
Random.seed!(SEED)

@testset "StochSystem" begin
    include("stochsystem.jl")
end

@testset "Large Deviations" begin
    include("largedeviations/MAM.jl")
    include("largedeviations/gMAM.jl")
end

@testset "Examples" begin
    include("examples/tutorial.jl")
end
