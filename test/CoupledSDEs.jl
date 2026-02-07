using CriticalTransitions, Test

@testset "API" begin
    using DynamicalSystemsBase, Test
    using OrdinaryDiffEq: Tsit5, OrdinaryDiffEqCore
    using StochasticDiffEq: SDEProblem, SRA, SOSRA, LambaEM, CorrelatedWienerProcess
    using CriticalTransitions

    @testset "drift" begin
        using LinearAlgebra
        using CriticalTransitions.CTLibrary: fitzhugh_nagumo

        σ = 0.2
        param = [1.0, 3, 1, 1, 1, 0]
        u0 = zeros(2)
        sys = CoupledSDEs(fitzhugh_nagumo, u0, param; noise_strength=σ)

        using CriticalTransitions: drift
        drift_vector = drift(sys, [0, 0])
        drift_vector isa SVector{2,Float64}
        @test drift_vector == [0, 0]
    end

    # Creation of lorenz
    @inbounds function lorenz_rule(u, p, t)
        σ = p[1]
        ρ = p[2]
        β = p[3]
        du1 = σ * (u[2] - u[1])
        du2 = u[1] * (ρ - u[3]) - u[2]
        du3 = u[1] * u[2] - β * u[3]
        return SVector{3}(du1, du2, du3)
    end
    @inbounds function lorenz_rule_iip(du, u, p, t)
        σ = p[1]
        ρ = p[2]
        β = p[3]
        du[1] = σ * (u[2] - u[1])
        du[2] = u[1] * (ρ - u[3]) - u[2]
        du[3] = u[1] * u[2] - β * u[3]
        return nothing
    end

    σ = 0.1
    function diagonal_noise!(σ)
        function (du, u, p, t)
            du .= σ .* ones(length(u))
            return nothing
        end
    end
    diagonal_noise(σ) = (u, p, t) -> SVector{3}(σ, σ, σ)

    u0 = [0, 10.0, 0]
    p0 = [10, 28, 8 / 3]
    Γ = [1.0 0.3 0.0; 0.3 1 0.5; 0.0 0.5 1.0]

    # diagonal additive noise
    lor_oop = CoupledSDEs(lorenz_rule, u0, p0)
    lor_iip = CoupledSDEs(
        SDEProblem(lorenz_rule_iip, diagonal_noise!(σ), copy(u0), (0.0, Inf), p0)
    )
    lor_SRA = CoupledSDEs(lorenz_rule, u0, p0; diffeq=(alg=SRA(), abstol=1e-2, reltol=1e-2))

    diffeq_cov = (alg=LambaEM(), abstol=1e-2, reltol=1e-2, dt=0.1)
    lor_oop_cov = CoupledSDEs(lorenz_rule, u0, p0; covariance=Γ, diffeq=diffeq_cov)
    lor_iip_cov = CoupledSDEs(lorenz_rule_iip, u0, p0; covariance=Γ, diffeq=diffeq_cov)

    @testset "correct SDE propagation" begin
        u0 = [0, 10.0, 0]
        p0 = [10, 28, 8 / 3]
        lorenz_oop = CoupledSDEs(lorenz_rule, u0, p0)
        @test lorenz_oop.integ.alg isa SOSRA

        lorenz_SRA = CoupledSDEs(
            lorenz_rule, u0, p0; diffeq=(alg=SRA(), abstol=1e-3, reltol=1e-3, verbose=false)
        )
        @test lorenz_SRA.integ.alg isa SRA

        # also test SDEproblem creation
        prob = lorenz_SRA.integ.sol.prob

        ds = CoupledSDEs(prob, (alg=SRA(), abstol=0.0, reltol=1e-3, verbose=false))

        @test ds.integ.alg isa SRA

        @test_throws ArgumentError CoupledSDEs(prob; diffeq=(alg=SRA(),))

        # CoupledODEs creation
        ds = CoupledODEs(lorenz_oop)
        @test dynamic_rule(ds).f == lorenz_rule
        @test ds.integ.alg isa Tsit5
        # and back
        sde = CoupledSDEs(ds, p0)
        @test dynamic_rule(sde).f.f == lorenz_rule
        @test sde.integ.alg isa SOSRA
    end

    @testset "interface" begin
        f(u, p, t) = 1.01u
        f!(du, u, p, t) = du .= 1.01u
        @testset "covariance" begin
            g(u, p, t) = sqrt([1 0.3; 0.3 1])
            corr = CoupledSDEs(f, zeros(2); covariance=[1 0.3; 0.3 1])
            corr_alt = CoupledSDEs(f, zeros(2); g=g, noise_prototype=zeros(2, 2))
            @test corr.noise_type == corr_alt.noise_type
            @test all(
                corr.integ.g(zeros(2), (), 0.0) .== corr_alt.integ.g(zeros(2), (), 0.0)
            )
        end

        @testset "ArgumentError" begin
            W = CorrelatedWienerProcess([1 0.3; 0.3 1], 0.0, zeros(2), zeros(2))
            @test_throws ArgumentError CoupledSDEs(f!, zeros(2); noise_process=W)

            g!(du, u, p, t) = du .= u
            @test_throws ArgumentError CoupledSDEs(
                f!, zeros(2); g=(g!), covariance=[1 0.3; 0.3 1]
            )

            g(u, p, t) = u
            @test_throws AssertionError CoupledSDEs(f!, zeros(2); g=g)

            Csde = CoupledSDEs(f!, zeros(2))
            diffeq = (alg=SRA(), abstol=1e-2, reltol=1e-2)
            @test_throws ArgumentError CoupledSDEs(Csde.integ.sol.prob; diffeq=diffeq)
        end
    end
end
