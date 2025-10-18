using CriticalTransitions
using DynamicalSystemsBase
using Test

@testset "RateSystem" begin
    function f(u, p, t) # out-of-place
        x = u[1]
        λ = p[1]
        dx = (x + λ)^2 - 1
        return SVector{1}(dx)
    end
    x0 = [-1.0]
    p_auto = [0.0]
    auto_sys = CoupledODEs(f, x0, p_auto)

    p(t) = tanh(t)         # A monotonic function that describes the parameter shift
    interval = (-5, 5)
    rc = RateConfig(p, interval)

    # 2) Then specify how you would like to use the RateConfig to ramp the pidx'th parameter of auto_sys:
    pidx = 1                # Index of the parameter-container of auto_sys that you want to ramp
    forcing_start = -50.0    # Time when the parameter shift should start (before this, the final system will be autonomous)
    forcing_length = 105.0   # Time-interval over which p(interval) is spread out or squeezed
    forcing_scale = 3.0     # Amplitude of the ramping. `p` is then automatically rescaled
    t0 = -70.0              # Initial time of the resulting non-autonomous system (relevant to later compute trajectories)

    rs = RateSystem(auto_sys, rc, pidx; forcing_start, forcing_length, forcing_scale, t0)

    @testset "obtaining info" begin
        @test current_state(rs) == x0
        @test DynamicalSystemsBase.initial_state(rs) == x0
        @test current_parameters(rs) == p_auto
        @test initial_parameters(rs) == p_auto
        @test current_time(rs) == t0
        @test initial_time(rs) == t0
        @test isinplace(rs) == false
        @test isdeterministic(rs) == false
        @test isdiscretetime(rs) == false
        @test rs(0) == x0
    end

    # Compute trajectories
    T = forcing_length + 40.0                      # length of the trajectory that we want to compute
    auto_traj = trajectory(auto_sys, T, x0)
    nonauto_traj = trajectory(rs.system, T, x0)

    @test isapprox(auto_traj[1][end, 1], -1; atol=1e-2)
    @test isapprox(nonauto_traj[1][end, 1], -4; atol=1e-2)
end

@testset "Rate system utils" begin
    rc = RateConfig(:linear)
    data = show(rc)
    @test data[1][1] == 0.0
    @test data[2][end] == 1.0
end
