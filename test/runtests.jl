using CriticalTransitions, StaticArrays
const CT = CriticalTransitions
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

@testset "Large Deviations" begin
    include("largedeviations/action_fhn.jl")
    include("largedeviations/MAM.jl")
    include("largedeviations/gMAM.jl")
end

# @testset "utilities" begin
#     include("utils.jl")
# end

# @testset "basin" begin
#     include("basin/basin_boundary.jl")
# end

@testset "Trajactories" begin
    include("trajactories/simulate.jl")
    include("trajactories/transition.jl")
end

@testset "Extentions" begin
    @testset "ChaosToolsExt" begin
        include("ext/ChaosToolsExt.jl")
    end

    @testset "CoupledSDEsBaisin" begin
        include("ext/CoupledSDEsBaisin.jl")
    end
end

@testset "Doctests" begin
    using Documenter
    Documenter.doctest(CriticalTransitions)
end
