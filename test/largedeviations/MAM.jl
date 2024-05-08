@testset "MAM FitzHugh-Nagumo" begin
    p = [0.1, 3, 1, 1, 1, 0]
    fhn = CoupledSDEs(fitzhugh_nagumo, diag_noise_funtion(0.1), zeros(2), p)
    x_i = SA[sqrt(2 / 3), sqrt(2 / 27)]
    x_f = SA[0.001, 0.0]
    N, T = 200, 2.0
    inst = min_action_method(
        fhn, x_i, x_f, N, T; maxiter = 100, verbose = false, save_info = false)
    S = fw_action(fhn, inst, range(0.0, T, length = N))
    @test isapprox(S, 0.18, atol = 0.01)
end

@testset "MAM Ornstein-Uhlenbeck" begin
    ou = CoupledSDEs((u, p, t) -> -u, diag_noise_funtion(0.18), SA[1.0])
    x0 = -1
    xT = 2.0
    T = 10.0
    N = 51
    t = range(0, T, N)
    inst_mam = min_action_method(
        ou, SA[x0], SA[xT], N, T, verbose = false, save_info = false)
    inst_sol = ((xT - x0 * exp(-T)) * exp.(t) .+ (x0 * exp(T) - xT) * exp.(-t)) /
               (exp(T) - exp(-T))
    @test maximum(abs.(inst_mam' .- inst_sol)) < 0.1
end
