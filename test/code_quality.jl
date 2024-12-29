@testset "Code quality" begin
    using ExplicitImports, Aqua
    ignore_deps = [:Random, :LinearAlgebra, :Printf, :Test, :Pkg]

    @test check_no_stale_explicit_imports(CriticalTransitions) == nothing
    @test check_all_explicit_imports_via_owners(CriticalTransitions) == nothing
    Aqua.test_ambiguities(CriticalTransitions)
    Aqua.test_all(
        CriticalTransitions;
        deps_compat=(
            ignore=ignore_deps,
            check_extras=(ignore=ignore_deps,),
            check_weakdeps=(ignore=ignore_deps,),
        ),
        piracies=(
            treat_as_own=[
                CriticalTransitions.DynamicalSystemsBase.SciMLBase.AbstractSDEIntegrator,
                CriticalTransitions.DynamicalSystemsBase.StateSpaceSets.StateSpaceSet
            ],
        ),
        ambiguities=false,
    )
end

@testset "Code linting" begin
    using JET
    JET.test_package(CriticalTransitions; target_defined_modules=true)
end
