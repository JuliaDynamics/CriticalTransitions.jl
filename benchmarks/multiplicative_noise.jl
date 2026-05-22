using CriticalTransitions, BenchmarkTools
using LinearAlgebra
using Random
# StaticArrays re-exported by CriticalTransitions; SA macro available via that.

const SUITE = BenchmarkGroup()

let
    function meier_stein(u, p, t)
        x, y = u
        return SA[x - x^3 - 10 * x * y^2, -(1 + x^2) * y]
    end
    ds = CoupledSDEs(meier_stein, zeros(2); noise_strength = 0.25)
    Nt = 100
    xx = range(-1.0, 1.0; length = Nt)
    yy = 0.3 .* (-xx .^ 2 .+ 1)
    x_initial = Matrix([xx yy]')
    SUITE["sgMAM/additive/meier_stein"] = @benchmarkable begin
        minimize_geometric_action(
            $(FreidlinWentzellHamiltonian(ds)), $x_initial,
            GeometricGradient(; stepsize = 1.0);
            maxiters = 100, show_progress = false,
        )
    end
end

let
    α = 0.3
    b1d(u, p, t) = SA[-u[1]]
    g1d(u, p, t) = SA[sqrt(1 + α * u[1]^2);;]
    ds = CoupledSDEs(b1d, SA[1.0]; g = g1d, noise_prototype = SMatrix{1, 1}(0.0))
    Nt = 80
    x_initial = reshape(collect(range(1.0, -1.0; length = Nt)), 1, Nt)
    SUITE["sgMAM/diagonal/1d_ou"] = @benchmarkable begin
        minimize_geometric_action(
            $(FreidlinWentzellHamiltonian(ds)), $x_initial,
            GeometricGradient(; stepsize = 1.0);
            maxiters = 100, show_progress = false,
        )
    end
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
    SUITE["sgMAM/general/2d_offdiag"] = @benchmarkable begin
        minimize_geometric_action(
            $(FreidlinWentzellHamiltonian(ds)), $x_initial,
            GeometricGradient(; stepsize = 1.0);
            maxiters = 100, show_progress = false,
        )
    end
end

SUITE
