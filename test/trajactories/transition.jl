
@testset "fitzhugh_nagumo" begin
    p = [1.0, 3.0, 1.0, 1.0, 1.0, 0.0] # Parameters (ϵ, β, α, γ, κ, I)
    σ = 0.24 # noise strength

    # StochSystem
    sys = CoupledSDEs(fitzhugh_nagumo, diag_noise_funtion(σ), zeros(2), p)

    # Calculate fixed points
    ds = CoupledODEs(sys)
    box = intervals_to_box([-2, -2], [2, 2])
    eqs, eigs, stab = fixedpoints(ds, box)

    # Store the two stable fixed points
    fp1, fp2 = eqs[stab]

    transition(sys, fp1, fp2)
end
