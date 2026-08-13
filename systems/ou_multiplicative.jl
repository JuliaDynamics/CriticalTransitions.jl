"""
    ou_multiplicative_1d(α)

Canonical 1D Ornstein-Uhlenbeck SDE with multiplicative diagonal noise
``\\sigma(x) = \\sqrt{1 + \\alpha x^2}``. Drift is ``-x``.

Used as a fixture by tests and examples that exercise state-dependent diagonal
diffusion through the Freidlin-Wentzell machinery.
"""
function ou_multiplicative_1d(α)
    b(u, p, t) = SA[-u[1]]
    g(u, p, t) = SA[sqrt(1 + α * u[1]^2);;]
    return CriticalTransitions.CoupledSDEs(
        b, SA[1.0]; g = g, noise_prototype = SMatrix{1, 1}(0.0),
    )
end

"""
    linear_offdiag_2d_sde()

Canonical 2D linear-drift SDE with full off-diagonal multiplicative noise. The
diffusion matrix has both diagonal entries `1 + 0.2 u_i` and off-diagonal entries
proportional to the opposite component.

Used as a fixture for `GeneralNoise` paths through the Freidlin-Wentzell machinery.
"""
function linear_offdiag_2d_sde()
    b(u, p, t) = SA[-u[1], -u[2]]
    function g(u, p, t)
        s11 = 1 + 0.2 * u[1]; s22 = 1 + 0.2 * u[2]
        s12 = 0.3 * u[2];     s21 = 0.3 * u[1]
        return SMatrix{2, 2}(s11, s21, s12, s22)
    end
    return CriticalTransitions.CoupledSDEs(
        b, SA[1.0, 0.0]; g = g, noise_prototype = SMatrix{2, 2}(zeros(2, 2)),
    )
end
