@testset "Simulate FitzHugh-Nagumo model" begin
    # System setup
    p = [0.1, 3.0, 1.0, 1.0, 1.0, 0.0] # Parameters (ϵ, β, α, γ, κ, I)
    σ = 0.1
    init = SA[1.0, 0.0]
    sys = CoupledSDEs(fitzhugh_nagumo, init, p; noise_strength=σ)
    traj = trajectory(sys, 10, init)
    relax = deterministic_orbit(sys, 10, init)
    
    @test traj[1][1,1] == 1.0
    @test isapprox(relax[1][end,1], 0.816; atol=1e-2)
    # These tests could be improved - Reyk
end
