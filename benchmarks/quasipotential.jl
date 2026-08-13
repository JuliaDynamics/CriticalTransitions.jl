function benchmark_quasipotential!(SUITE)
    # 2D non-gradient (additive noise): exercises the edge simplex pass and the
    # SMatrix `_line_integral` fast path.
    maier_stein(u, p, t) = SA[
        u[1] - u[1]^3 - 5 * u[1] * u[2]^2,
        -(1 + u[1]^2) * u[2],
    ]
    sys_ms = CoupledSDEs(maier_stein, [-1.0, 0.0]; noise_strength = 0.3)
    SUITE["Large deviation"]["Quasipotential"]["OLIM, Maier-Stein 121x81"] =
        @benchmarkable quasipotential(
        $sys_ms,
        CartesianGrid((-1.5, 1.5, 121), (-1.0, 1.0, 81)),
        [-1.0, 0.0]; show_progress = false,
    ) seconds = 30

    # 3D gradient well: exercises the triangle simplex pass with its FD Newton loop.
    quad3(u, p, t) = SA[-u[1], -u[2], -u[3]]
    sys_3d = CoupledSDEs(quad3, [0.0, 0.0, 0.0]; noise_strength = 1.0)
    SUITE["Large deviation"]["Quasipotential"]["OLIM, 3D quadratic 11x11x11"] =
        @benchmarkable quasipotential(
        $sys_3d,
        CartesianGrid((-1.0, 1.0, 11), (-1.0, 1.0, 11), (-1.0, 1.0, 11)),
        [0.0, 0.0, 0.0]; show_progress = false, band_radius = 2,
    ) seconds = 10

    # 2D multiplicative noise: exercises the `_QInvDynamic` callable and the
    # fallback `_line_integral` method (Qinv recomputed per Simpson node).
    b2(u, p, t) = SA[-u[1], -u[2]]
    function g2(u, p, t)
        s11 = 1 + 0.2 * u[1]; s22 = 1 + 0.2 * u[2]
        s12 = 0.3 * u[2];     s21 = 0.3 * u[1]
        return @SMatrix [s11 s12; s21 s22]
    end
    sys_mult = CoupledSDEs(b2, SA[0.0, 0.0]; g = g2, noise_prototype = SMatrix{2, 2}(zeros(2, 2)))
    SUITE["Large deviation"]["Quasipotential"]["OLIM, 2D multiplicative 81x81"] =
        @benchmarkable quasipotential(
        $sys_mult,
        CartesianGrid((-1.0, 1.0, 81), (-1.0, 1.0, 81)),
        [0.0, 0.0]; show_progress = false,
    ) seconds = 30
    return nothing
end
