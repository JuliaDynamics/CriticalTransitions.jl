
@testset "fitzhugh_nagumo" begin
    using Random, LinearAlgebra
    Random.seed!(SEED)

    p = [1.0, 3.0, 1.0, 1.0, 1.0, 0.0] # Parameters (ϵ, β, α, γ, κ, I)
    σ = 0.215 # noise strength

    # CoupledSDEs
    sys = CoupledSDEs(fitzhugh_nagumo, zeros(2), p; noise_strength=σ, seed=SEED)

    fp1 = [0.816, 0.272]
    fp2 = [-0.816, -0.272]
    trajectory, time, succes = CT.transition(sys, fp1, fp2)
    @test succes
    @test time[end] < 1e3
    @test norm(trajectory[:, 1] - fp1) < 0.1
    @test norm(trajectory[:, end] - fp2) < 0.1

    ensemble = transitions(sys, fp1, fp2, 10)
    @test ensemble.success_rate ≈ 1.0
    @test ensemble.t_trans ≈ 4.493941793363376 atol = 1e-2
    # SEED is different on github
    # SEED doesn;t work on github
    @test length(ensemble.times) == 11 broken = true
    @test ensemble.t_res ≈ 5299.98 broken = true
end
