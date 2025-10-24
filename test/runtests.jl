using CriticalTransitions, StaticArrays
const CT = CriticalTransitions
using Test

using Random
const SEED = 0xd8e5d8df
Random.seed!(SEED)

using CriticalTransitions.CTLibrary: fitzhugh_nagumo

"""
@testset "Code Quality" begin
    include("code_quality.jl")
end

@testset "CoupledSDEs" begin
    include("CoupledSDEs.jl")
    include("covariance.jl")
end

#@testset "ModelingToolkit" begin
#    include("ModelingToolkit.jl")
#end

@testset "Large Deviations" begin
    include("largedeviations/action_fhn.jl")
    include("largedeviations/MAM.jl")
    include("largedeviations/gMAM.jl")
    include("largedeviations/sgMAM.jl")
    include("largedeviations/string_method.jl")
    include("largedeviations/Maier_stein.jl")
    include("largedeviations/interpolate.jl")
end

@testset "Transition Path Theory" begin
    include("transition_path_theory.jl")
end

@testset "Utilities" begin
    include("utils.jl")
end

@testset "Trajectories" begin
    include("trajectories/simulate.jl")
    include("trajectories/transition.jl")
end
"""

@testset "R-tipping" begin
    #include("r_tipping/RateSystem.jl")
    include("r_tipping/RateSystem_debug.jl")
end

"""
@testset "Extensions" begin
    @testset "ChaosToolsExt" begin
        include("ext/ChaosToolsExt.jl")
    end

    @testset "DistMesh2D" begin
        include("ext/DistMesh2D.jl")
    end
end

@testset "Doctests" begin
    using Documenter
    Documenter.doctest(CriticalTransitions)
end
"""