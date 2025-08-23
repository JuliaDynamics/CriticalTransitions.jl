using CriticalTransitions
using Test

@testset "RateSystem" begin
    function f(u, p, t) # out-of-place
        x = u[1]
        dx = (x+p[1])^2 - 1
        return SVector{1}(dx)
    end

    x0 = [-1.0]
    auto_sys = CoupledODEs(f, x0, [0.0])

    function p(p_parameters, t)
        p_max = p_parameters[1]
        p_ = (p_max/2)*(tanh(p_max*t/2)+1)
        return SVector{1}(p_)
    end

    p_max = 3.0
    p_parameters = [p_max] # parameter of the function p

    r = 4/3-0.02   # r just below critical rate
    t_start = -Inf # start time of non-autonomous part
    t_end = Inf    # end time of non-autonomous part

    rp = CriticalTransitions.RateConfig(p, p_parameters, r, t_start,t_end)

    t0 = -10.0      # initial time of the system
    nonauto_sys = apply_ramping(auto_sys, rp, t0)

    T = 20.0        # final simulation time
    auto_traj = trajectory(auto_sys, T, x0)
    nonauto_traj = trajectory(nonauto_sys, T, x0)

    @test isapprox(auto_traj[1][end, 1], -1; atol=1e-4)
    @test isapprox(nonauto_traj[1][end, 1], -4; atol=1e-4)
end
