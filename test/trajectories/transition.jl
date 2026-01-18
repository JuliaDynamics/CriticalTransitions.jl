
@testset "fitzhugh_nagumo" begin
    using Random, LinearAlgebra
    Random.seed!(SEED)
    # SEED is different on github
    # SEED doesn;t work on github

    p = [1.0, 3.0, 1.0, 1.0, 1.0, 0.0] # Parameters (ϵ, β, α, γ, κ, I)
    σ = 0.315 # noise strength

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

    @test stats.success_rate == 1.0
    @test stats.transition_time < 10

    @test length(ensemble.times) == 10
    @test stats.residence_time < 100

    @testset "fitzhugh_nagumo - Nmax warning" begin
        # Test warning when Nmax is reached
        # Use very small tmax to make transitions fail
        p = [1.0, 3.0, 1.0, 1.0, 1.0, 0.0] # Parameters (ϵ, β, α, γ, κ, I)
        σ = 0.215 # noise strength

        # CoupledSDEs
        sys = CoupledSDEs(fitzhugh_nagumo, zeros(2), p; noise_strength=σ, seed=SEED)
        fp1 = [0.816, 0.272]
        fp2 = [-0.816, -0.272]

        @test_warn "Maximum number of attempts " ensemble_fail = transitions(
            sys, fp1, fp2, 10; Nmax=3, tmax=0.01, cut_start=false, show_progress=false
        )
    end
end
