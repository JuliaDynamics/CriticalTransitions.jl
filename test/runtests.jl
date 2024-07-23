using CriticalTransitions, StaticArrays
const CT = CriticalTransitions
using Test

using Random
const SEED = 0xd8e5d8df
Random.seed!(SEED)

function fitzhugh_nagumo(u, p, t)
    x, y = u
    ϵ, β, α, γ, κ, I = p

    dx = (-α * x^3 + γ * x - κ * y + I) / ϵ
    dy = -β * y + x

    return SA[dx, dy]
end

@testset "Code quality" begin
    using ExplicitImports, Aqua
    ignore_deps = [:Random, :LinearAlgebra, :Printf, :Test, :Pkg]

    @test check_no_stale_explicit_imports(CriticalTransitions) == nothing
    @test check_all_explicit_imports_via_owners(CriticalTransitions) == nothing
    # Aqua.test_ambiguities(CriticalTransitions)
    # Aqua.test_all(
    #     CriticalTransitions;
    #     deps_compat=(
    #         ignore=ignore_deps,
    #         check_extras=(ignore=ignore_deps,),
    #         check_weakdeps=(ignore=ignore_deps,),
    #     ),
    #     piracies=(treat_as_own=[],),
    #     ambiguities=false,
    # )
end

@testset "Code linting" begin
    using JET
    # JET.test_package(HarmonicBalance; target_defined_modules=true)
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
