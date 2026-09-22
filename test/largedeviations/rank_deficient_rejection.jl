using CriticalTransitions, StaticArrays
using Test

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

    sys = FreidlinWentzellHamiltonian(ds)
    Nt = 20
    xx = range(-1.0, 1.0; length = Nt)
    yy = 0.3 .* (-xx .^ 2 .+ 1)
    path = Matrix([xx yy]')

    # Both public geometric minimizers reject rank-deficient diffusion.
    err_s = try
        minimize_geometric_action(sys, path; maxiters = 1, show_progress = false); nothing
    catch e
        e
    end
    @test err_s isa ArgumentError
    @test occursin("rank-deficient", sprint(showerror, err_s))

    err_g = try
        minimize_geometric_action(ds, path; maxiters = 1, show_progress = false); nothing
    catch e
        e
    end
    @test err_g isa ArgumentError
    @test occursin("rank-deficient", sprint(showerror, err_g))
end

@testset "Static-vector diagonal diffusion rank check (#339)" begin
    drift(u, p, t) = SA[-u[1], -u[2]]
    diffusion(u, p, t) = SA[0.3 + 0.1 * u[1]^2, 0.4 + 0.1 * u[2]^2]
    ds = CoupledSDEs(drift, SA[0.1, -0.2]; g = diffusion)

    sys = FreidlinWentzellHamiltonian(ds)
    path = [range(-0.5, 0.5; length = 20)'; range(0.25, -0.25; length = 20)']

    res = minimize_geometric_action(
        sys, path, GeometricGradient(; stepsize = 1.0);
        maxiters = 1, show_progress = false,
    )
    @test res isa MinimumActionPath
    @test isfinite(res.action)
end
