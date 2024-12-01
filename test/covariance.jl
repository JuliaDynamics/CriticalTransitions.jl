StochasticSystemsBase = Base.get_extension(
    DynamicalSystemsBase, :StochasticSystemsBase
)
diffusion_matrix = StochasticSystemsBase.diffusion_matrix


@testset "diffusion_matrix" begin
    Γ = [1.0 0.3 0.0; 0.3 1 0.5; 0.0 0.5 1.0]
    A = sqrt(Γ)
    lorenz_oop = CoupledSDEs(lorenz_rule, u0, p0; covariance=Γ, diffeq=diffeq_cov)
    @test A ≈ diffusion_matrix(lorenz_oop)
    @test A isa AbstractMatrix
    @test Γ ≈ A * A'
end

@testset "approximate cov" begin
    @testset "StochasticDiffEq noise_prototype" begin
        f(u, p, t) = [0.0, 0.0]
        g(u, p, t) = sqrt([1.0 0.3; 0.3 1])
        noise_prototype = zeros(2, 2)

        sdeprob = SDEProblem(f, g, zeros(2), (0.0, 100_000), noise_rate_prototype=noise_prototype)
        sol = solve(sdeprob, EM(), saveat=0.1, dt=0.1)
        approx = cov(diff(reduce(hcat, sol.u); dims=2)./sqrt(0.1); dims=2)
        @test approx ≈ Γ atol = 1e-1
    end
    @testset "CT" begin
        using StatsBase: cov
        diffeq_cov = (alg=EM(), abstol=1e-2, reltol=1e-2, dt=0.1)

        Γ = [1.0 0.3; 0.3 1]
        f(u, p, t) = [0.0, 0.0]
        ds = CoupledSDEs(f, zeros(2), (); covariance=Γ, diffeq=diffeq_cov)

        tr, _ = trajectory(ds, 1_000, Δt=0.1)
        approx = cov(diff(reduce(hcat, tr.data); dims=2)./sqrt(0.1); dims=2)
        @test approx ≈ Γ atol = 1e-1 broken=true
    end
end

@testset "parametrize cov" begin
    Γ = [1.0 0.3; 0.3 1]
    f(u, p, t) = [0.0, 0.0]
    function diffusion(u, p, t)
        Γ = [1.0 p[1]; p[1] 1.0]
        return sqrt(Γ)
    end
    ds = CoupledSDEs(f, zeros(2), [0.3]; g=diffusion, noise_prototype=zeros(2, 2))
    A = diffusion_matrix(ds)
    @test Γ ≈ A * A'
    set_parameter!(ds, 1, 0.5)
    A = diffusion_matrix(ds)
    @test A * A' ≈ [1.0 0.5; 0.5 1.0]
end

@testset "noise" begin
    diff = diffusion_matrix(sys)
    cov = covariance_matrix(sys)
    @test diff isa AbstractMatrix
    @test cov isa AbstractMatrix
    @test all(diag(diff) .== σ)
    @test diff .^ 2 == cov

    W = noise_process(sys)
    int_cov = W.covariance
    # The internal covariance in the DiffEqNoiseProcess.NoiseProcess should be nothing
    @test int_cov == nothing
end
