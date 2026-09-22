using Random
Random.seed!(SEED)

p = [1.0, 3.0, 1.0, 1.0, 1.0, 0.0]
σ = 0.2
sys = CoupledSDEs(fitzhugh_nagumo, zeros(2), p; noise_strength = σ)

T, N = 2.0, 100
x_i = SA[sqrt(2 / 3), sqrt(2 / 27)]
x_f = SA[0.0, 0.0]
path = reduce(hcat, range(x_i, x_f; length = N))
time = range(0.0, T; length = N)

@testset "fw_action" begin
    S = fw_action(sys, path, time)
    @test isapprox(S, 0.32, atol = 0.01)
end

@testset "om_action" begin
    S = om_action(sys, path, time, σ)
    @test isapprox(S, 0.21, atol = 0.01)
end

@testset "action" begin
    @test action(sys, path, time, "FW") == fw_action(sys, path, time)
end

@testset "geometric_action" begin
    S = geometric_action(sys, path)
    @test isapprox(S, 0.23, atol = 0.01)
end

@testset "div_drift" begin
    @test div_drift(sys, zeros(2)) == -2.0
    @test div_drift(sys, x_i) == -4.0
end
