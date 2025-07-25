using CriticalTransitions
using Test

@testset "RateSystem" begin
    function f(u,p,t) # out-of-place
        x = u[1]
        λ = p[1]
        dx = (x+λ)^2 - 1
        return SVector{1}(dx)
    end

    lambda = 0.0 
    p = [lambda]
    x0 = [-1.]
    auto_sys = CoupledODEs(f,x0,p)

    function λ(p,t)
        λ_max = p[1]
        lambda = (λ_max/2)*(tanh(λ_max*t/2)+1)
        return SVector{1}(lambda)
    end

    λ_max = 3.
    p_lambda = [λ_max] # parameter of the function lambda

    r = 4/3-0.02   # r just below critical rate
    t_start = -Inf # start time of non-autonomous part
    t_end = Inf    # end time of non-autonomous part

    rp = CriticalTransitions.RateProtocol(λ,p_lambda,r,t_start,t_end)

    t0 = -10.      # initial time of the system
    nonauto_sys = RateSystem(auto_sys,rp,t0)

    T = 20.        # final simulation time
    auto_traj = trajectory(auto_sys,T,x0)
    nonauto_traj = trajectory(nonauto_sys,T,x0)

    @test isapprox(auto_traj[1][end,1], -1; atol=1e-4)
    @test isapprox(nonauto_traj[1][end,1], -4; atol=1e-4)
end