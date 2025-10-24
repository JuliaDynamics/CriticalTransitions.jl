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

    p(t) = tanh(t)
    rc = ForcingProfile(p, (-5.0, 5.0))

    pidx = 1
    forcing_start_time = 20.0
    forcing_length = 105.0
    forcing_scale = 1.5

    rs = RateSystem(ds, rc, pidx; forcing_start_time, forcing_length, forcing_scale)

    println("diffeq")
    println(rs.system.diffeq)
    println(ds.diffeq)

    println("Parameter as a function of time")
    println([get_forcing(rs, i) for i in range(0,140,length=15)])

    # Compute trajectories
    T = forcing_length + forcing_start_time + 40.0
    auto_traj = trajectory(ds, T, x0; Δt=10)
    nonauto_traj = trajectory(rs, T, x0; Δt=10)

    println("Unforced trajectory")
    println(auto_traj[1][:,1])

    println("---")
    println("Forced trajectory:")
    println("Time")
    println([Vector(nonauto_traj[2][:])[i] for i in 1:15])
    println("State")
    println([nonauto_traj[1][i,1] for i in 1:15])
    println("Parameter (this should be the same as above but it's not!)")
    println([get_forcing(rs, nonauto_traj[2][i]) for i in 1:15])
    println("Drift (going crazy)")
    println([dynamic_rule(rs)(nonauto_traj[1][i,1], get_forcing(rs, Vector(nonauto_traj[2][:])[i]), Vector(nonauto_traj[2][:])[i]) for i in 1:15])

    @test isapprox(auto_traj[1][end, 1], -1; atol=1e-2)
    @test isapprox(nonauto_traj[1][end, 1], -4; atol=1e-2)
end