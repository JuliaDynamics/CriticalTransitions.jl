
@testset "fitzhugh_nagumo" begin
    using Random
    # const SEED = 0xd8e5d8df
    Random.seed!(SEED)

    p = [1.0, 3.0, 1.0, 1.0, 1.0, 0.0] # Parameters (ϵ, β, α, γ, κ, I)
    σ = 0.215 # noise strength

    # CoupledSDEs
    sys = CoupledSDEs(fitzhugh_nagumo, diag_noise_funtion(σ), zeros(2), p, seed = SEED)

    # Calculate fixed points
    ds = CoupledODEs(sys)
    box = intervals_to_box([-2, -2], [2, 2])
    eqs, eigs, stab = fixedpoints(ds, box)

    # Store the two stable fixed points
    fp1, fp2 = eqs[stab]
    @test fp1 ≈ -fp2

    trajectory, time, succes = CT.transition(sys, fp1, fp2)
    @test succes
    @test time[end] < 1e3
    @test norm(trajectory[:, 1] - fp1) < 0.1
    @test norm(trajectory[:, end] - fp2) < 0.1

    ensemble = transitions(sys, fp1, fp2, 10)
    @test ensemble.success_rate ≈ 1.0
    @test ensemble.t_trans ≈ 4.493941793363376
    @test ensemble.t_res ≈ 5449.261866107592 skip=true # SEED is not working on github
    @test length(ensemble.times) ≈ 11 skip=true # SEED is not working on github
end
