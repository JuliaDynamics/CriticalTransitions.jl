using Test
using CriticalTransitions
using StaticArrays
using LinearAlgebra: tr

# Verifies the FW convention adopted by the package:
# `Q` is canonicalized to `tr(Q) = D`, so the action is independent of
# `noise_strength` and invariant under orthogonal changes of basis. See the
# docstring of `fw_action` for the underlying derivation.

@testset "Rotation invariance with non-diagonal Q" begin
    ou(u, p, t) = -u

    Q_diag = [1.0 0.0; 0.0 4.0]
    sys_diag = CoupledSDEs(ou, SA[0.0, 0.0]; covariance = Q_diag, noise_strength = 1.0)

    # 45-degree rotation
    R = [cos(π / 4) -sin(π / 4); sin(π / 4) cos(π / 4)]
    Q_rot = R * Q_diag * R'
    sys_rot = CoupledSDEs(ou, SA[0.0, 0.0]; covariance = Q_rot, noise_strength = 1.0)

    N = 201
    time_arr = range(0.0, 1.0; length = N)
    path_orig = zeros(2, N)
    for (i, t) in enumerate(time_arr)
        path_orig[:, i] = [1.0 + t, 0.5 + 0.3 * sin(2π * t)]
    end
    path_rot = R * path_orig

    # The FW action is a coordinate-independent geometric quantity, so the same
    # physical SDE/path in two bases must give the same number. With the previous
    # L1-of-entries normalization, these differed by tr(Q)/||Q||_1 across bases.
    @test isapprox(
        fw_action(sys_diag, path_orig, time_arr),
        fw_action(sys_rot, path_rot, time_arr);
        rtol = 1.0e-10,
    )
    @test isapprox(
        geometric_action(sys_diag, path_orig),
        geometric_action(sys_rot, path_rot);
        rtol = 1.0e-10,
    )
end

@testset "Action is independent of noise_strength" begin
    ou(u, p, t) = -u

    Q = [1.0 0.0; 0.0 4.0]
    sys_a = CoupledSDEs(ou, SA[0.0, 0.0]; covariance = Q, noise_strength = 0.3)
    sys_b = CoupledSDEs(ou, SA[0.0, 0.0]; covariance = Q, noise_strength = 2.5)

    N = 201
    time_arr = range(0.0, 1.0; length = N)
    path = zeros(2, N)
    for (i, t) in enumerate(time_arr)
        path[:, i] = [1.0 + t, 0.5 + 0.3 * sin(2π * t)]
    end

    @test isapprox(
        fw_action(sys_a, path, time_arr),
        fw_action(sys_b, path, time_arr);
        rtol = 1.0e-12,
    )
end

@testset "1D Ornstein-Uhlenbeck: action matches closed-form Φ_FW" begin
    # dx = -x dt + σ dW. Closed-form minimum action from 0 to x_T over time T:
    #   S_min(T) = x_T^2 / (1 - exp(-2T)) →  x_T^2  as T → ∞.
    # In 1D the trace-normalized covariance is always 1, so the returned action
    # is the σ-independent rate function `Φ_FW = x_T^2` (in the T → ∞ limit).
    σ = 0.5
    ou_1d(u, p, t) = SA[-u[1]]
    sys = CoupledSDEs(ou_1d, SA[0.0]; noise_strength = σ)

    x_T = 2.0
    T = 10.0
    N = 2001
    time_arr = range(0.0, T; length = N)
    # Instanton: φ(t) = x_T·sinh(t)/sinh(T) (solves ϕ̈ - ϕ = 0, BCs ϕ(0)=0, ϕ(T)=x_T)
    path = reshape([x_T * sinh(t) / sinh(T) for t in time_arr], 1, :)

    S_action = fw_action(sys, path, time_arr)
    S_expected = x_T^2 / (1 - exp(-2 * T))
    @test isapprox(S_action, S_expected; rtol = 1.0e-4)

    # σ-independence sanity check
    sys2 = CoupledSDEs(ou_1d, SA[0.0]; noise_strength = 2.0)
    @test isapprox(fw_action(sys2, path, time_arr), S_action; rtol = 1.0e-12)
end
