@testset "Simulate FitzHugh-Nagumo model" begin
    # System setup
    p = [0.1, 3.0, 1.0, 1.0, 1.0, 0.0] # Parameters (ϵ, β, α, γ, κ, I)
    σ = 0.1
    # Simulate
    sys = CoupledSDEs(fitzhugh_nagumo, idfunc, SA[1.0, 0.0], σ, p)
    traj = trajectory(sys, 10)
    sys = CoupledSDEs(fitzhugh_nagumo, idfunc, SA[1.0, 0.0], σ, p)
    sim = simulate(sys, 10)
    @test traj[1][1, 1] == 1.0
    @test sim.u[1][1] == 1.0
    # These tests could be improved - Reyk
end
