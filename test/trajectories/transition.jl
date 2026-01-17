
@testset "fitzhugh_nagumo" begin
    using Random, LinearAlgebra
    Random.seed!(SEED)

    p = [1.0, 3.0, 1.0, 1.0, 1.0, 0.0] # Parameters (ϵ, β, α, γ, κ, I)
    σ = 0.215 # noise strength

    # CoupledSDEs
    sys = CoupledSDEs(fitzhugh_nagumo, zeros(2), p; noise_strength=σ, seed=SEED)

    fp1 = [0.816, 0.272]
    fp2 = [-0.816, -0.272]
    tr, time, success = CT.transition(sys, fp1, fp2)
    @test success
    @test time[end] < 1e3
    @test norm(tr[1, :] - fp1) < 0.1
    @test norm(tr[end, :] - fp2) < 0.1

    ensemble = transitions(sys, fp1, fp2, 10)
    stats = ensemble.stats
    @show stats.success_rate
    @show stats.transition_time
    @show stats.residence_time
    @test isapprox(stats.success_rate, 0.833; atol=1e-2) ||
        isapprox(stats.success_rate, 0.909; atol=1e-2)
    @test isapprox(stats.transition_time, 5.213; atol=1e-2) ||
        isapprox(stats.transition_time, 5.6512; atol=1e-2)
    # SEED is different on github
    # SEED doesn;t work on github
    @test length(ensemble.times) == 10
    @test isapprox(stats.residence_time, 346.5424; atol=1e-2) ||
        isapprox(stats.residence_time, 177.70; atol=1e-2)

    # Test warning when Nmax is reached
    # Use very small tmax to make transitions fail
    @test_logs (:warn, r"Maximum number of attempts.*reached") begin
        ensemble_fail = transitions(sys, fp1, fp2, 5; Nmax=3, tmax=0.01, show_progress=false)
    end
end
