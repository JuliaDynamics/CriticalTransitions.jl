using CriticalTransitions
using DynamicalSystemsBase
using Test

pidx = 1
forcing_start_time = 20.0
forcing_duration = 100.0
forcing_scale = 1.5
T = forcing_duration + 40.0
x0 = [-1.0]
p_auto = [0.0]

# Autonomous drift
function f(u, p, t) # out-of-place
    x = u[1]
    λ = p[1]
    dx = (x + λ)^2 - 1
    return SVector{1}(dx)
end

# Hard-coded non-autonomous drift
function fexpl(u, p, t) # out-of-place
    x = u[1]
    λ_0, scale, shift, rate = p
    λ = λ_0 + scale*(tanh(rate*(t-shift)) + 1)
    dx = (x + λ)^2 - 1
    return SVector{1}(dx)
end

ds = CoupledODEs(f, x0, p_auto)

profile(t) = tanh(t)
fp = ForcingProfile(profile, (-5.0, 5.0))

rs = RateSystem(ds, fp, pidx; forcing_start_time, forcing_duration, forcing_scale)

@testset "frozen_system" begin
    frozen = frozen_system(rs, rs.forcing.t0)
    unforced = rs.forcing.unforced_rule
    u = [2.0];
    t = 10.0
    du_frozen = dynamic_rule(frozen)(u, current_parameters(frozen), t)
    du_orig = dynamic_rule(ds)(u, current_parameters(ds), t)
    du_unforced = unforced(u, [rs.forcing.p0], t)

    @test du_frozen == du_orig
    @test du_frozen == du_unforced
end

@testset "RateSystem" begin
    @testset "DynamicalSystems API" begin
        @test current_state(rs) == x0
        @test DynamicalSystemsBase.initial_state(rs) == x0
        @test current_parameters(rs) == p_auto
        @test initial_parameters(rs) == p_auto
        @test current_time(rs) == 0.0
        @test initial_time(rs) == 0.0
        @test isinplace(rs) == false
        @test isdeterministic(rs) == true
        @test isdiscretetime(rs) == false
        @test rs(0) == x0
    end

    # Compute trajectories
    auto_traj = trajectory(ds, T, x0)
    nonauto_traj = trajectory(rs, T, x0)

    @test isapprox(auto_traj[1][end, 1], -1; atol=1e-2)
    @test isapprox(nonauto_traj[1][end, 1], -4; atol=1e-2)

    @testset "Equivalence with hard-coded" begin
        forcing_start_time=0.0
        forcing_duration=100.0
        forcing_scale=1.0
        fp = ForcingProfile(profile, (-20.0, 20.0))

        ds = CoupledODEs(f, x0, p_auto)
        sys_constructed = RateSystem(
            ds, fp, pidx; forcing_start_time, forcing_duration, forcing_scale
        )

        p_hardcoded = [p_auto[1], forcing_scale, forcing_duration/2, 40/forcing_duration]
        sys_hardcoded = CoupledODEs(fexpl, x0, p_hardcoded; t0=0.0)

        tr_constructed, _ = trajectory(sys_constructed, T, x0)
        tr_hardcoded, _ = trajectory(sys_hardcoded, T, x0)

        @test all(tr_constructed .≈ tr_hardcoded)
    end
end
