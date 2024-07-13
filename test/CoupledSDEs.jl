@testset "initialisation" begin
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
        prob = SDEProblem(meier_stein, diagonal_noise(σ), zeros(2), (0.0, Inf))
        sde = CoupledSDEs(meier_stein, diagonal_noise(σ), zeros(2))

        @test sde.integ.sol.prob.f == prob.f
        @test sde.integ.sol.prob.g == prob.g
        @test sde.integ.sol.prob.u0 == prob.u0
        @test sde.integ.sol.prob.tspan == prob.tspan
        @test sde.integ.sol.prob.p == prob.p
        @test sde.integ.sol.prob.noise_rate_prototype == prob.noise_rate_prototype
        @test sde.integ.sol.prob.noise == prob.noise
        @test sde.integ.sol.prob.seed == prob.seed
    end
    # @show sde.integ.W.dist
    # @show sde.integ.sol.prob

    @testset "Scalar noise Wiener" begin
        # Scalar noise Wiener
        # a single random variable is applied to all dependent variables
        W = WienerProcess(0.0, 0.0, 0.0)
        prob = SDEProblem(
            meier_stein, diagonal_noise(σ), zeros(2), (0.0, Inf); noise=W
        )
        sde = CoupledSDEs(meier_stein, diagonal_noise(σ), zeros(2); noise=W)

        @test sde.integ.sol.prob.f == prob.f
        @test sde.integ.sol.prob.g == prob.g
        @test sde.integ.sol.prob.u0 == prob.u0
        @test sde.integ.sol.prob.tspan == prob.tspan
        @test sde.integ.sol.prob.p == prob.p
        @test sde.integ.sol.prob.noise_rate_prototype == prob.noise_rate_prototype
        @test sde.integ.sol.prob.noise == prob.noise
        @test sde.integ.sol.prob.seed == prob.seed
    end

    @testset "Multiplicative noise" begin
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

        @test sde.integ.sol.prob.f == prob.f
        @test sde.integ.sol.prob.g == prob.g
        @test sde.integ.sol.prob.u0 == prob.u0
        @test sde.integ.sol.prob.tspan == prob.tspan
        @test sde.integ.sol.prob.p == prob.p
        @test sde.integ.sol.prob.noise_rate_prototype == prob.noise_rate_prototype
        @test sde.integ.sol.prob.noise == prob.noise
        @test sde.integ.sol.prob.seed == prob.seed
    end
end
