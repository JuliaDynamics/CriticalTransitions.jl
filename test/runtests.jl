using CriticalTransitions, StaticArrays
const CT = CriticalTransitions
using Test

using Random
const SEED = 0xd8e5d8df
Random.seed!(SEED)

using CriticalTransitions.CTLibrary: fitzhugh_nagumo

@testset "code quality" begin
    include("code_quality.jl")
end

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

@testset "Utilities" begin
    include("utils.jl")
end

@testset "Trajectories" begin
    include("trajectories/simulate.jl")
    #   include("trajectories/transition.jl")
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
