@testset "Method error" begin
    @test_throws MethodError basins()
    @test_throws MethodError basboundary()
    @test_throws MethodError basinboundary()
end

# Test for intervals_to_box
@testset "intervals_to_box" begin
    using IntervalArithmetic, ChaosTools
    bmin = [-2, -1, 0]
    bmax = [2, 1, 1]
    expected = (interval([-2, 2]...), interval([-1, 1]...), interval([0, 1]...))
    @test CriticalTransitions.intervals_to_box(bmin, bmax).v.data == expected
end

@testset "fitzhugh_nagumo" begin
    using Attractors, ChaosTools
    # Define systems
    system = CoupledSDEs(
        fitzhugh_nagumo, [2.0, 0.1], [0.1, 3.0, 1.0, 1.0, 1.0, 0.0]; noise_strength=0.215
    )

    # Set up domain
    A, B, C = [0.0, 0.0], [1.0, 0.0], [0.0, 1.0]
    box = CriticalTransitions.intervals_to_box([-2.0, -2.0], [2.00, 2.0])
    step = 0.1

    boas = basins(system, A, B, C, box; bstep=[step, step])
    X, Y, attractors, h = boas
    boundary = basinboundary(boas)
    @test length(attractors) == 2
    @test length(unique(h)) == 3
    @test h[end ÷ 2 + 1, end ÷ 2 + 1] == -1
    @test boundary[:, end ÷ 2 + 1] == zeros(2)
end
