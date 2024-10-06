@testset "API" begin
    using LinearAlgebra
    function fitzhugh_nagumo(u, p, t)
        x, y = u
        ϵ, β, α, γ, κ, I = p

        dx = (-α * x^3 + γ * x - κ * y + I) / ϵ
        dy = -β * y + x

        return SA[dx, dy]
    end

    σ = 0.2
    param = [1.0, 3, 1, 1, 1, 0]
    u0 = zeros(2)
    sys = CoupledSDEs(fitzhugh_nagumo, u0, param; noise_strength=σ)

    @testset "correlation" begin end

    @testset "noise" begin
        diff = diffusion_matrix(sys)
        cov = covariance_matrix(sys)
        @test diff isa AbstractMatrix
        @test cov isa AbstractMatrix
        @test all(diag(diff) .== σ)
        @test diff .^ 2 == cov

        W = noise_process(sys)
        int_cov = W.covariance
        # The internal covariance in the DiffEqNoiseProcess.NoiseProcess should be nothing
        @test int_cov == nothing
    end

    @testset "drift" begin
        using CriticalTransitions: drift
        drift_vector = drift(sys, [0, 0])
        drift_vector isa SVector{2,Float64}
        @test drift_vector == [0, 0]
    end
end
