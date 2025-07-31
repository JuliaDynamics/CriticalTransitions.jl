
@testset "moving_sinks" begin
    using ChaosTools
    # Dynamical system
    function fhn(u,p,t)
        eps, b = p 
        x, y = u
        dx = (-x^3 + x - y)/eps
        dy = -b*y + x
        return SA[dx, dy]
    end

    p = [1,6]
    sys = CoupledODEs(fhn, [1.,1], p)

    # Forcing
    function λ(p, t)
        λ_max = p[2]
        lambda = (λ_max/2)*(tanh(λ_max*t/2)+1)
        return SVector{2}([p[1], lambda])
    end
    
    r = 0.1
    rp = CriticalTransitions.RateProtocol(λ, p, r, -10, 10)

    # Calculate moving sinks
    box = [interval(-2,2), interval(-1,1)]
    fp, eig, stab = moving_sinks(sys, rp, box; times=0:0.1:1)
    @test length(fp[1]) == 3
end