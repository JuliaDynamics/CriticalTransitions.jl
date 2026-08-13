using CriticalTransitions
using Attractors
using Test
using Random

# 1D double well with an additive forced parameter: ẋ = x - x³ + λ.
# At λ = 0 it is bistable (x ≈ ±1); ramping λ up destroys the left well, so the
# track/return/tip outcomes of Ritchie2023 all appear for a range of rates.
dwell(u, p, t) = SVector(u[1] - u[1]^3 + p[1])

function rate_diagram(t0)
    ds = CoupledODEs(dwell, [-1.0], [0.0]; t0 = t0)
    grid = (range(-2.0, 2.0; length = 101),)
    mapper = AttractorsViaRecurrences(ds, grid; sparse = true)
    sampler, = statespace_sampler(grid)
    profile = ForcingProfile(x -> cos(x)^2, (-π / 2, 0.0))   # 0 -> 1 ramp
    rs = RateSystem(ds, profile, 1; forcing_start_time = 20.0, reverse = true, t0 = t0)
    Δps = range(0.0, 0.9; length = 6)
    Δts = 2.0 .^ range(-1, 3; length = 6)
    Random.seed!(1)
    rate_type, attractors_cont = rate_track_return_tip(
        rs, Δts, Δps, mapper, sampler; proximity_kw = (ε = 0.1,), u0 = [-1.0]
    )
    return rate_type, attractors_cont
end

@testset "rate_track_return_tip" begin
    rate_type, attractors_cont = rate_diagram(0.0)

    @test rate_type isa Matrix{Int}
    @test size(rate_type) == (6, 6)                  # length(Δps) × length(Δts)
    @test length(attractors_cont) == 6               # one entry per Δp
    # tracking, return-but-not-track, and tipping all occur.
    @test sort(unique(rate_type)) == [1, 2, 3]

    # `forcing_start_time` is absolute and `reinit!` resets to t0, so the stepping in
    # `rate_track_return_tip` must reach the forcing-interval boundaries independently
    # of t0. The whole diagram should be invariant under a shift of t0.
    rate_type_shifted, _ = rate_diagram(5.0)
    @test rate_type_shifted == rate_type
end
