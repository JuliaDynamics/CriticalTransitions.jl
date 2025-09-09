using CriticalTransitions
using Test

@testset "RateSystem" begin
    function f(u,p,t) # out-of-place
        x = u[1]
        λ = p[1]
        dx = (x+λ)^2 - 1
        return SVector{1}(dx)
    end;
    x0 = [-1.0];
    p_auto = [0.];  
    auto_sys = CoupledODEs(f,x0,p_auto);
    

    pidx=1;                 # Index of the parameter that is made nonautonomous 
    p(t) = tanh(t);         # Function that describes the parameter shift
    section_start = -100;   # start time of the section of the tanh we want to consider
    section_end = 100;      # start time of the section of the tanh we want to consider
    window_start = -50.;    # time when the parameter shift should start (before this, the system is autonomous)
    window_length = 105;    # tracks!
    dp=3;                   # amplitude of the parameter shift

    rc = CriticalTransitions.RateConfig(pidx,p,section_start, section_end, window_start, window_length, dp);
    t0 = window_start - 20.0;      # initial time of the system
    nonauto_sys = apply_ramping(auto_sys, rc, t0);
    
    T = window_length + 40.0;        # simulation time
    auto_traj = trajectory(auto_sys, T, x0)   
    nonauto_traj = trajectory(nonauto_sys, T, x0)

    @test isapprox(auto_traj[1][end, 1], -1; atol=1e-2)
    @test isapprox(nonauto_traj[1][end, 1], -4; atol=1e-2)
end
