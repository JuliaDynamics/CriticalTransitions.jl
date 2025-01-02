@testset "autonomous invertible additive noise" begin
    f!(du, u, p, t) = du .= 1.01u
    g_corr_alt!(du, u, p, t) = (du .= [1 0.3 1; 0.3 1 1]; return nothing)
    corr_alt = CoupledSDEs(f!, zeros(2); g = g_corr_alt!, noise_prototype = zeros(2, 3))

    g_addit_non_autom!(du, u, p, t) = (du .= 1 / (1 + t); return nothing)
    addit_non_autom = CoupledSDEs(f!, zeros(2); g = g_addit_non_autom!)

    g_linear_multipli!(du, u, p, t) = (du .= u; return nothing)
    linear_multipli = CoupledSDEs(f!, rand(2) ./ 10; g = g_linear_multipli!)

    x_i = SA[sqrt(2 / 3), sqrt(2 / 27)]
    x_f = SA[0.001, 0.0]
    N, T = 200, 2.0
    init = reduce(hcat, range(x_i, x_f; length=N))

    for inst in [corr_alt, addit_non_autom, linear_multipli]
        @test_throws ArgumentError min_action_method(inst,init, T)
        @test_throws ArgumentError geometric_min_action_method(inst, init)
        @test_throws ArgumentError sgmam(inst, init)
    end
end
