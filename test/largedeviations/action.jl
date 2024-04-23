using Random
Random.seed!(SEED)

p = [1.0, 3.0, 1.0, 1.0, 1.0, 0.0] # Parameters (ϵ, β, α, γ, κ, I)
σ = 0.215 # noise strength

# StochSystem
sys = CoupledSDEs(fitzhugh_nagumo, diag_noise_funtion(σ), zeros(2), p, seed = SEED)

# Calculate fixed points
ds = CoupledODEs(sys)
box = intervals_to_box([-2, -2], [2, 2])
eqs, eigs, stab = fixedpoints(ds, box)

# Store the two stable fixed points
fp1, fp2 = eqs[stab]
@test fp1 ≈ -fp2

trajectory, time, succes = CT.transition(sys, fp1, fp2)
CT.fw_action(sys, trajectory, time)
trajectory
A = inv(CT.covariance_matrix(sys))

integrand = fw_integrand(sys, path, time, A)

# Test fw_action function
@testset "fw_action" begin

    path = rand(2, 10)
    time = collect(1:10)
    expected = 0.0 # Replace with expected value
    @test fw_action(sys, path, time) ≈ expected
end

# Test om_action function
@testset "om_action" begin
    sys = CoupledSDEs()
    path = rand(2, 10)
    time = collect(1:10)
    expected = 0.0 # Replace with expected value
    @test om_action(sys, path, time) ≈ expected
end

# Test action function
@testset "action" begin
    sys = CoupledSDEs()
    path = rand(2, 10)
    time = collect(1:10)
    expected_fw = 0.0 # Replace with expected value for fw_action
    expected_om = 0.0 # Replace with expected value for om_action
    @test action(sys, path, time, "FW") ≈ expected_fw
    @test action(sys, path, time, "OM") ≈ expected_om
end

# Test geometric_action function
@testset "geometric_action" begin
    sys = CoupledSDEs()
    path = rand(2, 10)
    arclength = 1.0
    expected = 0.0 # Replace with expected value
    @test geometric_action(sys, path, arclength) ≈ expected
end

# Test fw_integrand function
@testset "fw_integrand" begin
    sys = CoupledSDEs()
    path = rand(2, 10)
    time = collect(1:10)
    A = rand(2, 2)
    expected = zeros(10) # Replace with expected values
    @test fw_integrand(sys, path, time, A) ≈ expected
end

# Test div_b function
@testset "div_b" begin
    sys = CoupledSDEs()
    x = rand(2)
    expected = 0.0 # Replace with expected value
    @test div_b(sys, x) ≈ expected
end

# Test path_velocity function
@testset "path_velocity" begin
    path = rand(2, 10)
    time = collect(1:10)
    expected = zeros(2, 10) # Replace with expected values
    @test path_velocity(path, time) ≈ expected
end
