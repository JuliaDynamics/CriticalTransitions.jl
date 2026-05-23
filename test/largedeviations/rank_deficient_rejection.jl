using CriticalTransitions, StaticArrays
using Test

const CT = CriticalTransitions

@testset "Rank-deficient rejection (#325 deferred)" begin
    function langevin(u, p, t)
        x, p_ = u
        return SA[p_, -x - 0.1 * p_]
    end
    g_langevin(u, p, t) = SA[0.0 0.0; 0.0 sqrt(0.2)]
    ds = CoupledSDEs(
        langevin, SA[0.0, 0.0]; g = g_langevin,
        noise_prototype = SMatrix{2, 2}(zeros(2, 2)),
    )

    # Constructor accepts the Hamiltonian (no a(x) sampling at construction).
    sys = FreidlinWentzellHamiltonian(ds)
    @test sys isa FreidlinWentzellHamiltonian

    # Rejection happens at cache build (sgMAM) or workspace build (gMAM).
    Nt = 20
    xx = range(-1.0, 1.0; length = Nt)
    yy = 0.3 .* (-xx .^ 2 .+ 1)
    path = Matrix([xx yy]')

    err_s = try
        CT.build_sgmam_cache(sys, path, Nt); nothing
    catch e
        e
    end
    @test err_s isa ArgumentError
    @test occursin("rank-deficient", err_s.msg)
    @test occursin("FreidlinWentzellHamiltonian", err_s.msg)

    err_g = try
        minimize_geometric_action(ds, path); nothing
    catch e
        e
    end
    @test err_g isa ArgumentError
    @test occursin("rank-deficient", err_g.msg)
end
