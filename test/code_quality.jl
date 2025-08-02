using CriticalTransitions, Test

@testset "best practices" begin
    using Aqua

    Aqua.test_ambiguities(CriticalTransitions)
    Aqua.test_all(
        CriticalTransitions;
        piracies=(
            treat_as_own=[
                CriticalTransitions.DynamicalSystemsBase.SciMLBase.AbstractSDEIntegrator,
                CriticalTransitions.DynamicalSystemsBase.StateSpaceSets.StateSpaceSet,
            ],
        ),
        ambiguities=false,
    )
end

@testset "ExplicitImports" begin
    using ExplicitImports

    @test check_no_implicit_imports(CriticalTransitions) == nothing
    @test check_all_explicit_imports_via_owners(CriticalTransitions) == nothing
    @test check_all_explicit_imports_are_public(CriticalTransitions) == nothing
    @test check_no_stale_explicit_imports(CriticalTransitions) == nothing
    @test check_all_qualified_accesses_via_owners(CriticalTransitions) == nothing
    @test isnothing(
        check_all_qualified_accesses_are_public(
            CriticalTransitions;
            skip=(Base => Base.Experimental, Base => Core),
            ignore=(
            :EnsembleAlgorithm,
            :Terminated,
            :covariance_matrix,
            :diffusion_matrix,
            :get_extension,
            :register_error_hint,
        ),
        ),
    )
    @test check_no_self_qualified_accesses(CriticalTransitions) == nothing
end

if isempty(VERSION.prerelease)
    @testset "Code linting" begin
        using JET
        JET.test_package(CriticalTransitions; target_defined_modules=true)
    end
end
