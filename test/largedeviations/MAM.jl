using CriticalTransitions
using Test

using CriticalTransitions.CTLibrary: fitzhugh_nagumo

@testset "MAM FitzHugh-Nagumo" begin
    p = [0.1, 3, 1, 1, 1, 0]
    σ = 0.1
    fhn = CoupledSDEs(fitzhugh_nagumo, zeros(2), p; noise_strength=σ)
    x_i = SA[sqrt(2 / 3), sqrt(2 / 27)]
    x_f = SA[0.001, 0.0]
    N, T = 200, 2.0
    inst = min_action_method(fhn, x_i, x_f, T; N, maxiter=500, show_progress=false)
    # If you evolve for longer the path splits into two :/
    S = fw_action(fhn, Matrix(inst.path)', range(0.0, T; length=N))
    @test isapprox(S, 0.18, atol=0.01)
end

@testset "MAM Ornstein-Uhlenbeck" begin
    σ = 0.18
    ou = CoupledSDEs((u, p, t) -> -u, SA[1.0]; noise_strength=σ)
    x0 = -1
    xT = 2.0
    T = 10.0
    N = 51
    t = range(0, T, N)
    inst_mam = min_action_method(
        ou, SA[x0], SA[xT], T; N=51, maxiter=2000, show_progress=false
    )
    inst_sol =
        ((xT - x0 * exp(-T)) * exp.(t) .+ (x0 * exp(T) - xT) * exp.(-t)) /
        (exp(T) - exp(-T))
    @test maximum(abs.(inst_mam.path[:, 1] .- inst_sol)) < 0.1
end
