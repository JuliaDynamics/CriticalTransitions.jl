using Test

@testset "Interpolations" begin
    using CriticalTransitions: interpolate_path!
    using LinearAlgebra
    @testset "1D path" begin
        path = collect(reshape(range(0, 1, 10), 1, 10))
        α = zeros(10)
        s = range(0, 1; length=10)

        interpolate_path!(path, α, s)

        @test vec(path) ≈ α atol = 1e-10
        @test α ≈ s atol = 1e-10
    end

    @testset "2D path" begin
        path = [0.0 1.0 0.0; 0.0 1.0 2.0]
        α = zeros(3)
        s = range(0, 1; length=3)

        original = copy(path)
        interpolate_path!(path, α, s)

        # Check endpoints remain unchanged
        @test path[:, 1] ≈ original[:, 1]
        @test path[:, end] ≈ original[:, end]

        # Check if points are equally spaced
        diffs = vec(sqrt.(sum(diff(path; dims=2) .^ 2; dims=1)))
        @test all(x -> isapprox(x, diffs[1]; rtol=1e-10), diffs)

        # Check if α represents normalized cumulative distances
        @test α ≈ [0.0, 0.5, 1.0] atol = 1e-10
    end
end
