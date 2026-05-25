using Test
using CriticalTransitions
using PeriodicOrbits
using StaticArrays

@testset "LimitCycleFrame from PeriodicOrbit" begin
    function sl(u, p, t)
        x, y = u
        μ, ω = p
        r2 = x^2 + y^2
        return SA[(μ - r2) * x - ω * y, (μ - r2) * y + ω * x]
    end
    sys = CoupledSDEs(sl, [1.0, 0.0], (1.0, 1.0); noise_strength = 0.1)
    Nτ = 100
    T = 2π
    pts = StateSpaceSet([SA[cos(t), sin(t)] for t in range(0, T; length = Nτ + 1)[1:Nτ]])
    po = PeriodicOrbit(pts, T, true)
    lc = LimitCycleFrame(po, sys)
    @test lc.period ≈ T
    @test size(lc.E) == (2, 2, Nτ)
end
