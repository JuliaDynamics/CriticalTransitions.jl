using CriticalTransitions
using Test

using LinearAlgebra

const CT = CriticalTransitions

@testset "MultipleShooting rejects CoupledSDEs" begin
    ou1d(u, p, t) = SA[-u[1]]
    ds = CoupledSDEs(ou1d, [0.0]; noise_strength = 1.0)
    x_init = Matrix(reshape(collect(range(0.0, 1.0; length = 10)), 1, 10))
    err = try
        minimize_geometric_action(ds, x_init, MultipleShooting())
        nothing
    catch e
        e
    end
    @test err isa ArgumentError
    msg = sprint(showerror, err)
    @test occursin("FreidlinWentzellHamiltonian", msg)
    @test occursin("shooting", lowercase(msg))
end

@testset "Free endpoint rejected" begin
    bistable(u, p, t) = SA[u[1] - u[1]^3]
    ds = CoupledSDEs(bistable, [0.0]; noise_strength = 1.0)
    H = FreidlinWentzellHamiltonian(ds)
    # Start at a free (non-fixed) point.
    x_init = Matrix(reshape(collect(range(0.5, 1.0; length = 10)), 1, 10))
    err = try
        minimize_geometric_action(H, x_init, MultipleShooting(; nshoots = 4))
        nothing
    catch e
        e
    end
    @test err isa ArgumentError
    msg = sprint(showerror, err)
    @test occursin("fixed-point endpoint", msg)
    @test occursin("GeometricGradient", msg)
end

@testset "Non-hyperbolic fixed point rejected" begin
    # x' = -x^3 + y, y' = -x. Linearization at origin: J = [0 1; -1 0], purely imaginary eigenvalues.
    function rotrhs(u, p, t)
        x, y = u
        return SA[-x^3 + y, -x]
    end
    ds = CoupledSDEs(rotrhs, [0.0, 0.0]; noise_strength = 1.0)
    H = FreidlinWentzellHamiltonian(ds)
    x_init = Matrix(zeros(2, 10))
    err = try
        minimize_geometric_action(H, x_init, MultipleShooting(; nshoots = 4))
        nothing
    catch e
        e
    end
    @test err isa ArgumentError
    @test occursin("hyperbolic", sprint(showerror, err))
end
