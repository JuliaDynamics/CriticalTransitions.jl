
using ChaosTools

function meier_stein(u, p, t) # out-of-place
    x, y = u
    dx = x - x^3 - 10 * x * y^2
    dy = -(1 + x^2) * y
    return SA[dx, dy]
end
σ = 0.1

sde = CoupledSDEs(meier_stein, zeros(2), (); noise_strength=σ)

fps, eigs, stab = fixedpoints(sde, [-3, -3], [3, 3])

@test stab == [true, true, false]
fp1, fp2 = fps[stab]
@test fp1 ≈ -fp2
@test all(broadcast(v -> all(v .< 0), real.(eigs[stab])))
