@testset "initialisation" begin
    function diagonal_noise!(σ)
        function (du, u, p, t)
            idfunc!(du, u, p, t)
            du .*= σ
            return nothing
        end
    end
    diagonal_noise(σ) = (u, p, t) -> σ .* idfunc(u, p, t)

    # Define the system
    function meier_stein(u, p, t) # out-of-place
        x, y = u
        dx = x - x^3 - 10 * x * y^2
        dy = -(1 + x^2) * y
        return SA[dx, dy]
    end
    σ = 0.1

    @testset "diagonal additive noise" begin
        # diagonal additive noise: σ*N(0,dt)
        # a vector of random numbers dW whose size matches the output of g where the noise is applied element-wise
        prob = SDEProblem(meier_stein, diagonal_noise(σ), zeros(2), (0.0, Inf), ())
        sde = CoupledSDEs(meier_stein, idfunc, zeros(2), (), σ)

        @test sde.integ.sol.prob.f isa SDEFunction
        @test sde.integ.sol.prob.f.f.f == prob.f.f
        # @test sde.integ.sol.prob.g == prob.g
        @test sde.integ.sol.prob.u0 == prob.u0
        @test sde.integ.sol.prob.tspan == prob.tspan
        @test sde.integ.sol.prob.p == prob.p
        @test sde.integ.sol.prob.noise_rate_prototype == prob.noise_rate_prototype
        @test sde.integ.sol.prob.noise == prob.noise
        @test sde.integ.sol.prob.seed == prob.seed

        noise = sde.integ.sol.prob.noise
        @test noise == nothing
    end
    # @show sde.integ.W.dist
    # @show sde.integ.sol.prob

    @testset "Scalar noise Wiener" begin
        # Scalar noise Wiener
        # a single random variable is applied to all dependent variables
        W = WienerProcess(0.0, 0.0, 0.0)
        prob = SDEProblem(meier_stein, diagonal_noise(σ), zeros(2), (0.0, Inf), (); noise=W)
        sde = CoupledSDEs(meier_stein, idfunc, zeros(2), (), σ; noise=W)

        @test sde.integ.sol.prob.f isa SDEFunction
        @test sde.integ.sol.prob.f.f.f == prob.f.f
        # @test sde.integ.sol.prob.g == prob.g
        @test sde.integ.sol.prob.u0 == prob.u0
        @test sde.integ.sol.prob.tspan == prob.tspan
        @test sde.integ.sol.prob.p == prob.p
        @test sde.integ.sol.prob.noise_rate_prototype == prob.noise_rate_prototype
        @test sde.integ.sol.prob.noise == prob.noise
        @test sde.integ.sol.prob.seed == prob.seed

        noise = sde.integ.sol.prob.noise
        @test noise == W
        @test length(noise.dW) == 1
        @test W.covariance == nothing
    end

    @testset "multiplicative noise Wiener" begin
        # multiplicative noise Wiener
        # a single random variable is applied to all dependent variables
        g_sde(u, p, t) = σ .* u
        g_Csde(u, p, t) = σ .* u
        prob = SDEProblem(meier_stein, g_sde, zeros(2), (0.0, Inf), ())
        sde = CoupledSDEs(meier_stein, g_Csde, zeros(2), (), σ)

        @test sde.integ.sol.prob.f isa SDEFunction
        @test sde.integ.sol.prob.f.f.f == prob.f.f
        @test sde.integ.sol.prob.u0 == prob.u0
        @test sde.integ.sol.prob.tspan == prob.tspan
        @test sde.integ.sol.prob.p == prob.p
        @test sde.integ.sol.prob.noise_rate_prototype == prob.noise_rate_prototype
        @test sde.integ.sol.prob.noise == prob.noise
        @test sde.integ.sol.prob.seed == prob.seed

        noise = sde.integ.sol.prob.noise
        @test noise == nothing
    end

    @testset "Non-diagonal noise" begin
        # non-diagonal noise allows for the terms to linearly mixed via g being a matrix.
        # four Wiener processes and two dependent random variables
        # the output of g to be a 2x4 matrix, such that the solution is g(u,p,t)*dW, the matrix multiplication.
        f(du, u, p, t) = du .= 1.01u
        function g(du, u, p, t)
            du[1, 1] = 0.3u[1]
            du[1, 2] = 0.6u[1]
            du[1, 3] = 0.9u[1]
            du[1, 4] = 0.12u[1]
            du[2, 1] = 1.2u[2]
            du[2, 2] = 0.2u[2]
            du[2, 3] = 0.3u[2]
            du[2, 4] = 1.8u[2]
            return nothing
        end
        # prob = SDEProblem(f, g, ones(2), (0.0, 1.0), noise_rate_prototype = zeros(2, 4))
        # diffeq =(alg=SRIW1(), noise_rate_prototype = zeros(2, 4))
        sde = CoupledSDEs(
            f, g, zeros(2); noise_rate_prototype=zeros(2, 4), diffeq=(alg=RKMilCommute(),)
        )
        prob = SDEProblem(f, g, zeros(2), (0.0, Inf); noise_rate_prototype=zeros(2, 4))

        @test sde.integ.sol.prob.f isa SDEFunction
        @test sde.integ.sol.prob.f.f.f == prob.f.f
        # @test sde.integ.sol.prob.g == prob.g
        @test sde.integ.sol.prob.u0 == prob.u0
        @test sde.integ.sol.prob.tspan == prob.tspan
        @test sde.integ.sol.prob.p == prob.p
        @test sde.integ.sol.prob.noise_rate_prototype == prob.noise_rate_prototype
        @test sde.integ.sol.prob.noise == prob.noise
        @test sde.integ.sol.prob.seed == prob.seed
    end

    @testset "Correlated noise Wiener" begin
        f!(du, u, p, t) = du .= 1.01u

        ρ = 0.3
        Γ = [1 ρ; ρ 1]
        t0 = 0.0
        W0 = zeros(2)
        Z0 = zeros(2)
        W = CorrelatedWienerProcess(Γ, t0, W0, Z0)
        prob = SDEProblem(f!, diagonal_noise!(σ), zeros(2), (0.0, Inf), (); noise=W)
        sde = CoupledSDEs(f!, idfunc!, zeros(2), (), σ; noise=W)

        @test sde.integ.sol.prob.f isa SDEFunction
        @test sde.integ.sol.prob.f.f.f == prob.f.f
        @test sde.integ.sol.prob.u0 == prob.u0
        @test sde.integ.sol.prob.tspan == prob.tspan
        @test sde.integ.sol.prob.p == prob.p
        @test sde.integ.sol.prob.noise_rate_prototype == prob.noise_rate_prototype
        @test sde.integ.sol.prob.noise == prob.noise
        @test sde.integ.sol.prob.seed == prob.seed

        @test W.covariance == Γ
        @testset "covariance" begin
            using DiffEqNoiseProcess, StatsBase
            prob = NoiseProblem(W, (0.0, 1.0))
            output_func = (sol, i) -> (sol.dW, false)
            ensemble_prob = EnsembleProblem(prob; output_func=output_func)

            dt = 0.1
            sol = Array(solve(ensemble_prob; dt=dt, trajectories=1_000_000))

            @test zero(mean(sol; dims=2)[:]) ≈ mean(sol; dims=2)[:] atol = 1e-2
            @test W.covariance ≈ cov(sol; dims=2) / dt rtol = 1e-2
        end
    end
end
