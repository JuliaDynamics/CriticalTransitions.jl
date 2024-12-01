# Test for intervals_to_box
@testset "intervals_to_box" begin
    using IntervalArithmetic
    bmin = [-2, -1, 0]
    bmax = [2, 1, 1]
    expected = (interval([-2, 2]...), interval([-1, 1]...), interval([0, 1]...))
    @test intervals_to_box(bmin, bmax).v.data == expected
end

# Test for anorm
@testset "anorm" begin
    using CriticalTransitions: anorm
    vector = [1, 2, 3]
    A = [1 0 0; 0 2 0; 0 0 3]
    expected = sqrt(1^2 + 2^2 * 2 + 3^2 * 3)
    @test anorm(vector, A) == expected
end

# Test for subnorm
@testset "subnorm" begin
    using CriticalTransitions: subnorm
    vector = [3, 7, 4]
    directions = [1, 3]
    expected = sqrt(3^2 + 4^2)
    @test subnorm(vector; directions=directions) == expected
end