using ChaosTools

@testset "intervals_to_box" begin
    using ChaosTools
    using IntervalArithmetic: isequal_interval
    bmin = [-2, -1, 0]
    bmax = [2, 1, 1]
    expected = [interval([-2, 2]...), interval([-1, 1]...), interval([0, 1]...)]

    @test all(isequal_interval.(CriticalTransitions.intervals_to_box(bmin, bmax), expected))
end

@testset "fixedpoints Maier-Stein model" begin
    function maier_stein(u, p, t)
        x, y = u
        dx = x - x^3 - 10 * x * y^2
        dy = -(1 + x^2) * y
        return SA[dx, dy]
    end
    σ = 0.1

    sde = CoupledSDEs(maier_stein, zeros(2), (); noise_strength=σ)

    fps, eigs, stab = fixedpoints(sde, [-3, -3], [3, 3])

    @test stab == [true, false, true]
    fp1, fp2 = fps[stab]
    @test fp1 ≈ -fp2
    @test all(broadcast(v -> all(v .< 0), real.(eigs[stab])))
end