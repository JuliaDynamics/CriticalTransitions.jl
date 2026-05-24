@testset "autonomous invertible additive noise" begin
    f!(du, u, p, t) = du .= 1.01u
    g_corr_alt!(du, u, p, t) = (du .= [1 0.3 1; 0.3 1 1]; return nothing)
    corr_alt = CoupledSDEs(f!, zeros(2); g = (g_corr_alt!), noise_prototype = zeros(2, 3))

    g_addit_non_autom!(du, u, p, t) = (du .= 1 / (1 + t); return nothing)
    addit_non_autom = CoupledSDEs(f!, zeros(2); g = (g_addit_non_autom!))

    g_linear_multipli!(du, u, p, t) = (du .= u; return nothing)
    linear_multipli = CoupledSDEs(f!, rand(2) ./ 10; g = (g_linear_multipli!))

    x_i = SA[sqrt(2 / 3), sqrt(2 / 27)]
    x_f = SA[0.001, 0.0]
    N, T = 200, 2.0
    init = reduce(hcat, range(x_i, x_f; length = N))

    # `minimize_action` (MAM) and `minimize_geometric_action` (gMAM) both only
    # reject non-autonomous noise at the precheck stage; multiplicative /
    # state-dependent / correlated noise is supported (rank-deficient `a(x)` is
    # rejected later, at workspace / cache build).
    @test_throws ArgumentError minimize_action(addit_non_autom, init, T)
    @test_throws ArgumentError minimize_geometric_action(addit_non_autom, init)
end
