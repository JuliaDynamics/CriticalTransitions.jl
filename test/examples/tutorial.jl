@testset "Tutorial" begin
    function fitzhugh_nagumo(u, p, t)
        u, v = u
        ϵ, β, α, γ, κ, I = p[1]

        du = (-α * u^3 + γ * u - κ * v + I) / ϵ
        dv = -β * v + u

        SA[du, dv]
    end

    p = [1.0, 3.0, 1.0, 1.0, 1.0, 0.0] # Parameters (ϵ, β, α, γ, κ, I)
    σ = 0.18 # noise strength

    # StochSystem
    sys = StochSystem(fitzhugh_nagumo, p, zeros(2), σ)

    # Calculate fixed points
    ds = CoupledODEs(sys)
    box = intervals_to_box([-2, -2], [2, 2])
    eqs, eigs, stab = fixedpoints(ds, box)

    # Store the two stable fixed points
    fp1, fp2 = eqs[stab]

    sim = simulate(sys, fp1, dt = 0.01, tmax = 1e3)
end # @testset "Tutorial"
