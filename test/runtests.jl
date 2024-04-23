using CriticalTransitions, StaticArrays
CT = CriticalTransitions
using Test

using Random
const SEED = 0xd8e5d8df
Random.seed!(SEED)

@testset "CoupledSDEs" begin
    include("CoupledSDEs.jl")
end

@testset "ModelingToolkit" begin
    include("ModelingToolkit.jl")
end

# @testset "StochSystem" begin
#     include("stochsystem.jl")
# end

# @testset "Large Deviations" begin
#     include("largedeviations/MAM.jl")
#     include("largedeviations/gMAM.jl")
# end

# @testset "utilities" begin
#     include("utils.jl")
# end

# @testset "Baisin" begin
#     include("baisin/baisin_boundary.jl")
# end

@testset "Trajactories" begin
    include("trajactories/simulate.jl")
    include("trajactories/transition.jl")
end
