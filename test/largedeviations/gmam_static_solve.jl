using LinearAlgebra
using StaticArrays
using Test

@testset "static direct-gMAM covariance solve matches dense" begin
    a_static = @SMatrix [2.0 0.25; 0.25 1.5]
    b = [0.3, -0.7]
    φp = [1.2, 0.4]

    θ_static = zeros(2)
    Ainv_b_static = zeros(2)
    Ainv_φp_static = zeros(2)
    λ_static = CT._lambda_theta!(
        θ_static, a_static, b, φp, Ainv_b_static, Ainv_φp_static, nothing,
    )

    θ_dense = zeros(2)
    Ainv_b_dense = zeros(2)
    Ainv_φp_dense = zeros(2)
    λ_dense = CT._lambda_theta!(
        θ_dense, Matrix(a_static), b, φp, Ainv_b_dense, Ainv_φp_dense, nothing,
    )

    @test λ_static ≈ λ_dense
    @test θ_static ≈ θ_dense
    @test Ainv_b_static ≈ Ainv_b_dense
    @test Ainv_φp_static ≈ Ainv_φp_dense
end

@testset "static diagonal direct-gMAM covariance solve matches Diagonal" begin
    a_static = @SMatrix [2.0 0.0; 0.0 3.0]
    a_diag = Diagonal([2.0, 3.0])
    b = [0.4, -0.2]
    φp = [0.6, 0.8]

    θ_static = zeros(2)
    Ainv_b_static = zeros(2)
    Ainv_φp_static = zeros(2)
    λ_static = CT._lambda_theta!(
        θ_static, a_static, b, φp, Ainv_b_static, Ainv_φp_static, nothing,
    )

    θ_diag = zeros(2)
    Ainv_b_diag = zeros(2)
    Ainv_φp_diag = zeros(2)
    λ_diag = CT._lambda_theta!(
        θ_diag, a_diag, b, φp, Ainv_b_diag, Ainv_φp_diag, nothing,
    )

    @test λ_static ≈ λ_diag
    @test θ_static ≈ θ_diag
    @test Ainv_b_static ≈ Ainv_b_diag
    @test Ainv_φp_static ≈ Ainv_φp_diag
end
