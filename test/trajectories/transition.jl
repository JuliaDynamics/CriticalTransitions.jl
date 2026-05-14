@testset "fitzhugh_nagumo" begin
    using Random, LinearAlgebra
    Random.seed!(SEED)
    # SEED is different on github
    # SEED doesn;t work on github

    p = [1.0, 3.0, 1.0, 1.0, 1.0, 0.0] # Parameters (ϵ, β, α, γ, κ, I)
    σ = 0.315 # noise strength

    # CoupledSDEs
    sys = CoupledSDEs(fitzhugh_nagumo, zeros(2), p; noise_strength = σ, seed = SEED)

    fp1 = [0.816, 0.272]
    fp2 = [-0.816, -0.272]
    tr, time, success = CT.transition(sys, fp1, fp2)
    @test success
    @test time[end] < 1.0e3
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
        sys = CoupledSDEs(fitzhugh_nagumo, zeros(2), p; noise_strength = σ, seed = SEED)
        fp1 = [0.816, 0.272]
        fp2 = [-0.816, -0.272]

        @test_warn "Maximum number of attempts " ensemble_fail = transitions(
            sys, fp1, fp2, 10; Nmax = 3, tmax = 0.01, cut_start = false, show_progress = false
        )
    end

    @testset "sys.diffeq forwarded (#283)" begin
        # If maxiters from sys.diffeq is forwarded, the solver gives up before
        # reaching fp2 and all attempts fail.
        using StochasticDiffEq: SOSRA
        sys_low = CoupledSDEs(
            fitzhugh_nagumo, zeros(2), p;
            noise_strength = σ, seed = SEED,
            diffeq = (alg = SOSRA(), maxiters = 50),
        )
        ens = @test_logs (:warn, r"Maximum number of attempts") match_mode = :any transitions(
            sys_low, fp1, fp2, 3; Nmax = 3, show_progress = false
        )
        @test ens.stats.success_rate < 1
    end

    @testset "seed independence (#284)" begin
        sys_noseed = CoupledSDEs(fitzhugh_nagumo, zeros(2), p; noise_strength = σ)

        # Intra-call: no two parallel sims should share an RNG seed (the original bug).
        ens = transitions(sys_noseed, fp1, fp2, 3; show_progress = false)
        @test allunique(ens.times)

        # Cross-call: independent calls produce independent ensembles by default.
        e1 = transitions(sys_noseed, fp1, fp2, 2; show_progress = false)
        e2 = transitions(sys_noseed, fp1, fp2, 2; show_progress = false)
        @test e1.times != e2.times
    end

    @testset "seed kwarg reproducibility" begin
        sys_noseed = CoupledSDEs(fitzhugh_nagumo, zeros(2), p; noise_strength = σ)
        e1 = transitions(sys_noseed, fp1, fp2, 2; seed = 42, show_progress = false)
        e2 = transitions(sys_noseed, fp1, fp2, 2; seed = 42, show_progress = false)
        @test e1.times == e2.times

        tr1, t1, _ = CT.transition(sys_noseed, fp1, fp2; seed = 42)
        tr2, t2, _ = CT.transition(sys_noseed, fp1, fp2; seed = 42)
        @test t1 == t2
        @test tr1 == tr2
    end
end
