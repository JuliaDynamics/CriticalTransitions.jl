using LinearAlgebra
using StaticArrays
using Test

@testset "static coupled sgMAM covariance factor matches dense solve" begin
    a = @SMatrix [2.0 0.25; 0.25 1.5]
    rhs = [0.3, -0.7]
    expected = Matrix(a) \ rhs

    factor = CT._StaticCoupledCovarianceFactor(inv(a))
    work = copy(rhs)
    @test LinearAlgebra.ldiv!(factor, work) === work
    @test work ≈ expected
end

@testset "dynamic coupled sgMAM covariance keeps LU fallback" begin
    a = [2.0 0.25; 0.25 1.5]
    factor = CT._state_dependent_coupled_factor(a, copy(a))
    @test factor isa LinearAlgebra.LU

    rhs = [0.3, -0.7]
    expected = a \ rhs
    work = copy(rhs)
    LinearAlgebra.ldiv!(factor, work)
    @test work ≈ expected
end
