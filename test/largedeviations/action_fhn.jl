using Random
Random.seed!(SEED)

# System setup - FitzHugh-Nagumo model
p = [1.0, 3.0, 1.0, 1.0, 1.0, 0.0] # Parameters (ϵ, β, α, γ, κ, I)
σ = 0.2 # noise strength
sys = CoupledSDEs(fitzhugh_nagumo, zeros(2), p; noise_strength=σ)

A = inv(CriticalTransitions.normalize_covariance!(covariance_matrix(sys)))
T, N = 2.0, 100

x_i = SA[sqrt(2 / 3), sqrt(2 / 27)]
x_f = SA[0.0, 0.0]

path = reduce(hcat, range(x_i, x_f; length=N))
time = range(0.0, T; length=N)

@testset "fw_action" begin
    S = fw_action(sys, path, time)
    @test isapprox(S, 0.32, atol=0.01)
end

# Test om_action function
@testset "om_action" begin
    S = om_action(sys, path, time, σ)
    @test isapprox(S, 0.26, atol=0.01)
end

# Test action function
@testset "action" begin
    @test action(sys, path, time, "FW") == fw_action(sys, path, time)
end

# Test geometric_action function
@testset "geometric_action" begin
    S = geometric_action(sys, path)
    @test isapprox(S, 0.23, atol=0.01)
end

# Test fw_integrand function
@testset "fw_integrand" begin
    integrand = CriticalTransitions.fw_integrand(sys, path, time, A)
    @test all(integrand .>= -eps())
    # Integral of integrand should equal 2 * fw_action for same path/time.
    dt = time[2] - time[1]
    trapz = sum((integrand[2:end] .+ integrand[1:(end - 1)]) .* dt ./ 2)
    @test isapprox(trapz, 2 * fw_action(sys, path, time); rtol=1e-10, atol=1e-12)
end

# Test div_drift function
@testset "div_drift" begin
    @test CT.div_drift(sys, zeros(2)) == -2.0
    @test CT.div_drift(sys, x_i) == -4.0
end
