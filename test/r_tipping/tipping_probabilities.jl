@testset "Tipping probabilities" begin
    basins_before = [1 1 2; 1 2 2]
    basins_after = [1 2 2; -1 2 1]
    expected = [1 / 3 1 / 3 1 / 3; 1 / 3 2 / 3 0]

    P = tipping_probabilities(basins_before, basins_after)
    @test P == expected
    @test vec(sum(P; dims = 2)) == ones(2)

    # Labels need not be consecutive. Rows/columns follow sorted labels,
    # with the divergent label -1 placed last.
    noncontiguous_before = [5, 5, 2, 2]
    noncontiguous_after = [9, -1, 9, 3]
    @test tipping_probabilities(noncontiguous_before, noncontiguous_after) ==
        [1 / 2 1 / 2 0; 0 1 / 2 1 / 2]

    @test size(tipping_probabilities(Int[], Int[])) == (0, 0)
    @test_throws DimensionMismatch tipping_probabilities([1, 2], reshape([1, 2], 1, 2))

    grid = (range(0, 1; length = 2), range(0, 1; length = 3))
    attractors = Dict(k => StateSpaceSet([[Float64(k), 0.0]]) for k in 1:2)
    BoA_before = Attractors.ArrayBasinsOfAttraction(basins_before, attractors, grid)
    BoA_after = Attractors.ArrayBasinsOfAttraction(basins_after, attractors, grid)
    @test tipping_probabilities(BoA_before, BoA_after) == expected
end
