function benchmark_multiplicative_noise!(SUITE)
    let
        α = 0.3
        b1d(u, p, t) = SA[-u[1]]
        g1d(u, p, t) = SA[sqrt(1 + α * u[1]^2);;]
        ds = CoupledSDEs(b1d, SA[1.0]; g = g1d, noise_prototype = SMatrix{1, 1}(0.0))
        Nt = 80
        x_initial = reshape(collect(range(1.0, -1.0; length = Nt)), 1, Nt)
        sys = FreidlinWentzellHamiltonian(ds)
        opt = GeometricGradient(; stepsize = 1.0)
        SUITE["Large deviation"]["Simple geometric minimal action"]["multiplicative diagonal"] = @benchmarkable minimize_geometric_action(
            $sys, $x_initial, $opt; maxiters = 100, show_progress = false
        )
        SUITE["Large deviation"]["Geometric minimal action"]["multiplicative diagonal"] = @benchmarkable minimize_geometric_action(
            $ds, $x_initial, GeometricGradient(); maxiters = 100, show_progress = false
        ) seconds = 10
    end

    let
        b2(u, p, t) = SA[-u[1], -u[2]]
        function g2(u, p, t)
            s11 = 1 + 0.2 * u[1]; s22 = 1 + 0.2 * u[2]
            s12 = 0.3 * u[2];     s21 = 0.3 * u[1]
            return @SMatrix [s11 s12; s21 s22]
        end
        ds = CoupledSDEs(b2, SA[1.0, 0.0]; g = g2, noise_prototype = SMatrix{2, 2}(zeros(2, 2)))
        Nt = 60
        xx = collect(range(1.0, 0.0; length = Nt))
        yy = collect(range(0.0, 1.0; length = Nt))
        x_initial = Matrix([xx yy]')
        sys = FreidlinWentzellHamiltonian(ds)
        opt = GeometricGradient(; stepsize = 1.0)
        SUITE["Large deviation"]["Simple geometric minimal action"]["multiplicative off-diagonal"] = @benchmarkable minimize_geometric_action(
            $sys, $x_initial, $opt; maxiters = 100, show_progress = false
        )
        SUITE["Large deviation"]["Geometric minimal action"]["multiplicative off-diagonal"] = @benchmarkable minimize_geometric_action(
            $ds, $x_initial, GeometricGradient(); maxiters = 100, show_progress = false
        ) seconds = 10
    end
    return nothing
end
