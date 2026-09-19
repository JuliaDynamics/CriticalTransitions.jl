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
    time = collect(range(0.0, 1.0; length = Nt))
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

    SUITE["Large deviation"]["Action kernels"]["Onsager-Machlup"] = @benchmarkable om_action(
        $ds, $path, $time, 0.25
    ) seconds = 5

    α_noise = 0.3
    b_diag(u, p, t) = SA[-u[1]]
    g_diag(u, p, t) = SA[sqrt(1 + α_noise * u[1]^2);;]
    ds_diag = CoupledSDEs(
        b_diag, SA[1.0]; g = g_diag, noise_prototype = SMatrix{1, 1}(0.0),
    )
    path_diag = reshape(collect(range(1.0, -1.0; length = 80)), 1, 80)
    ws_diag = CT.geometric_gradient_workspace(ds_diag, path_diag)
    SUITE["Large deviation"]["Direct gMAM kernels"]["multiplicative diagonal step"] = @benchmarkable CT.geometric_gradient_step!(
        $ws_diag, $ds_diag, $path_diag; stepsize = 0.1
    ) seconds = 5

    b_coupled(u, p, t) = SA[-u[1], -u[2]]
    function g_coupled(u, p, t)
        s11 = 1 + 0.2 * u[1]
        s22 = 1 + 0.2 * u[2]
        s12 = 0.3 * u[2]
        s21 = 0.3 * u[1]
        return @SMatrix [s11 s12; s21 s22]
    end
    ds_coupled = CoupledSDEs(
        b_coupled, SA[1.0, 0.0];
        g = g_coupled, noise_prototype = SMatrix{2, 2}(zeros(2, 2)),
    )
    xx_coupled = collect(range(1.0, 0.0; length = 60))
    yy_coupled = collect(range(0.0, 1.0; length = 60))
    path_coupled = Matrix([xx_coupled yy_coupled]')
    ws_coupled = CT.geometric_gradient_workspace(ds_coupled, path_coupled)
    SUITE["Large deviation"]["Direct gMAM kernels"]["multiplicative off-diagonal step"] = @benchmarkable CT.geometric_gradient_step!(
        $ws_coupled, $ds_coupled, $path_coupled; stepsize = 0.1
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
