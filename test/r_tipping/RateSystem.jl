using CriticalTransitions
using Test

@testset "RateSystem" begin
    function f(u, p, t) # out-of-place
        x = u[1]
        λ = p[1]
        dx = (x+λ)^2 - 1
        return SVector{1}(dx)
    end;
    x0 = [-1.0];
    p_auto = [0.0];
    auto_sys = CoupledODEs(f, x0, p_auto);

    p(t) = tanh(t);         # A monotonic function that describes the parameter shift
    interval = (-5,5)
    rc = RateConfig(p, interval)

    # 2) Then specify how you would like to use the RateConfig to ramp the pidx'th parameter of auto_sys:
    pidx = 1                # Index of the parameter-container of auto_sys that you want to ramp
    forcing_start = -50.0    # Time when the parameter shift should start (before this, the final system will be autonomous)
    forcing_length = 105.0   # Time-interval over which p(interval) is spread out or squeezed
    forcing_scale = 3.0     # Amplitude of the ramping. `p` is then automatically rescaled 
    t0 = -70.0              # Initial time of the resulting non-autonomous system (relevant to later compute trajectories)

    RateSys = RateSystem(
        auto_sys,
        rc,
        pidx;
        forcing_start=forcing_start,
        forcing_length=forcing_length,
        forcing_scale=forcing_scale,
        t0=t0,
    )

    # Compute trajectories
    T = forcing_length + 40.0;                      # length of the trajectory that we want to compute
    auto_traj = trajectory(auto_sys, T, x0);
    nonauto_traj = trajectory(RateSys.system, T, x0);

    @test isapprox(auto_traj[1][end, 1], -1; atol=1e-2)
    @test isapprox(nonauto_traj[1][end, 1], -4; atol=1e-2)
end
