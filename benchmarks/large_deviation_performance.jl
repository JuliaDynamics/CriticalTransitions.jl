import CriticalTransitions as CT

function benchmark_large_deviation_performance!(SUITE)
    function maier_stein(u, p, t)
        x, y = u
        return SA[x - x^3 - 10 * x * y^2, -(1 + x^2) * y]
    end

    ds = CoupledSDEs(maier_stein, zeros(2); noise_strength = 0.25)
    fw = FreidlinWentzellHamiltonian(ds)

    Nt = 200
    xx = collect(range(-1.0, 1.0; length = Nt))
    yy = @. 0.3 * (1 - xx^2)
    path = Matrix([xx yy]')
    momentum = similar(path)
    @inbounds for j in axes(momentum, 2)
        momentum[1, j] = 0.05 * sinpi((j - 1) / (Nt - 1))
        momentum[2, j] = 0.03 * cospi((j - 1) / (Nt - 1))
    end
    Hp_buf = similar(path)
    Hx_buf = similar(path)

    SUITE["Large deviation"]["Hamiltonian kernels"]["auto H_p"] = @benchmarkable CT._eval_Hp!(
        $Hp_buf, $fw, $path, $momentum
    ) seconds = 5

    SUITE["Large deviation"]["Hamiltonian kernels"]["auto H_x"] = @benchmarkable CT._eval_Hx!(
        $Hx_buf, $fw, $path, $momentum
    ) seconds = 5

    init = path[:, 1:4:end]
    sg_opt = GeometricGradient(; max_backtracks = 0, stepsize = 1.0)
    SUITE["Large deviation"]["Fixed iteration"]["auto sgMAM 10 iterations"] = @benchmarkable minimize_geometric_action(
        $fw, $init, $sg_opt; maxiters = 10, show_progress = false
    ) seconds = 10

    gg_opt = GeometricGradient(; max_backtracks = 0, stepsize = 0.1)
    SUITE["Large deviation"]["Fixed iteration"]["direct gMAM 100 iterations"] = @benchmarkable minimize_geometric_action(
        $ds, $init, $gg_opt; maxiters = 100, show_progress = false
    ) seconds = 10

    return nothing
end
