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
