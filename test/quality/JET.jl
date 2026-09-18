using CriticalTransitions
using Test
using JET

@static if isempty(VERSION.prerelease)
    @testset "Code linting" begin
        JET.test_package(CriticalTransitions; target_modules = (CriticalTransitions,))
    end
end
