using LinearAlgebra
using StaticArrays
using Test

@testset "block Thomas solve matches dense block-tridiagonal reference" begin
    eps_step = 0.2
    diagonal_blocks = [
        [2.0 0.1; 0.1 1.6],
        [1.8 -0.15; -0.15 2.2],
        [2.4 0.2; 0.2 1.9],
    ]
    rhs = [0.3 1.2 -0.4; -0.7 0.4 0.9]
    rhs_reference = copy(rhs)
    # Include a zero interior lambda to certify that the unscaled recurrence
    # remains valid without dividing by lambda^2 or falling back.
    lambda = reshape([0.0, 1.1, 0.0, 1.3, 0.0], 1, :)

    Nx = 2
    L = length(diagonal_blocks)
    M = zeros(Nx * L, Nx * L)
    for i in 1:L
        rows = ((i - 1) * Nx + 1):(i * Nx)
        M[rows, rows] .= diagonal_blocks[i]
        if i < L
            next_rows = (i * Nx + 1):((i + 1) * Nx)
            q_i = lambda[i + 1]^2
            q_next = lambda[i + 2]^2
            for k in 1:Nx
                M[rows[k], next_rows[k]] = -eps_step * q_i
                M[next_rows[k], rows[k]] = -eps_step * q_next
            end
        end
    end
    expected = M \ vec(rhs_reference)

    schur = deepcopy(diagonal_blocks)
    factor_seeds = [Matrix{Float64}(I, Nx, Nx) for _ in 1:L]
    factors = [
        cholesky!(Hermitian(S, :L); check = false) for S in factor_seeds
    ]
    inv_prev = zeros(Nx, Nx)
    tmp = zeros(Nx)

    @test CT._block_thomas_solve!(
        schur, factors, rhs, inv_prev, tmp, eps_step, lambda,
    )
    @test vec(rhs) ≈ expected
end

@testset "coupled sgMAM selects block cache" begin
    drift(u, p, t) = SA[-u[1], -u[2]]
    function diffusion(u, p, t)
        return @SMatrix [1 + 0.1 * u[1] 0.2 * u[2]; -0.1 * u[1] 1 - 0.1 * u[2]]
    end
    ds = CoupledSDEs(
        drift, SA[0.2, -0.1];
        g = diffusion, noise_prototype = SMatrix{2, 2}(zeros(2, 2)),
    )
    sys = FreidlinWentzellHamiltonian(ds)
    path = [range(-0.4, 0.4; length = 12)'; range(0.2, -0.2; length = 12)']
    cache = CT.build_sgmam_cache(sys, path, size(path, 2))
    @test cache isa CT.SgMAMBlockCoupledCache
end
