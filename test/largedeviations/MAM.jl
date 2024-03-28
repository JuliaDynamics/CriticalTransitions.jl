@testset "OuMAM" begin
    ou = StochSystem((u, p, t) -> -u, [], [1.0], 1.0)
    x0 = -1
    xT = 2.0
    T = 10.0
    N = 51
    t = range(0, T, N)
    inst_mam = min_action_method(ou, [x0], [xT], N, T, verbose = false, save_info = false)
    inst_sol = ((xT - x0 * exp(-T)) * exp.(t) .+ (x0 * exp(T) - xT) * exp.(-t)) /
               (exp(T) - exp(-T))
    @test maximum(abs.(inst_mam' .- inst_sol)) < 0.1
end
