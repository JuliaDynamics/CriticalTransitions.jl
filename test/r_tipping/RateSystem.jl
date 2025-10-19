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
    ds = CoupledODEs(f, x0, p_auto)

    p(t) = tanh(t)         # A monotonic function that describes the parameter shift
    ptspan = (-5, 5)
    rc = RateFunction(p, ptspan)

    pidx = 1
    forcing_start_time = -50.0
    forcing_length = 105.0
    forcing_scale = 3.0
    t0 = -70.0

    rs = RateSystem(ds, rc, pidx; forcing_start_time, forcing_length, forcing_scale, t0)

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
    T = forcing_length + 40.0
    auto_traj = trajectory(ds, T, x0)
    nonauto_traj = trajectory(rs, T, x0)

    @test isapprox(auto_traj[1][end, 1], -1; atol=1e-2)
    @test isapprox(nonauto_traj[1][end, 1], -4; atol=1e-2)
end
