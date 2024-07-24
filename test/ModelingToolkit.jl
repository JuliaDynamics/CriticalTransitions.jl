using ModelingToolkit

@testset begin
    @independent_variables t
    D = Differential(t)
    sts = @variables x(t) y(t) z(t)
    ps = @parameters σ ρ
    @brownian β η
    s = 0.001
    β *= s
    η *= s

    eqs = [
        D(x) ~ σ * (y - x) + x * β,
        D(y) ~ x * (ρ - z) - y + y * β + x * η,
        D(z) ~ x * y - β * z + (x * z) * β,
    ]
    @named sys1 = System(eqs, t)
    sys1 = structural_simplify(sys1)

    drift_eqs = [D(x) ~ σ * (y - x), D(y) ~ x * (ρ - z) - y, D(z) ~ x * y]

    diffusion_eqs = [
         s*x 0
        s*y s*x
        (s * x * z)-s * z 0
    ]

    sys2 = SDESystem(drift_eqs, diffusion_eqs, t, sts, ps; name=:sys1)
    @test sys1 == sys2

    prob = SDEProblem(sys1, sts .=> [1.0, 0.0, 0.0], (0.0, 100.0), ps .=> (10.0, 26.0))
    solve(prob, LambaEulerHeun(); seed=1)
end

@testset begin
    using ModelingToolkit
    @independent_variables t
    ps = @parameters α β
    sts = @variables x(t)
    D = Differential(t)

    eqs = [D(x) ~ α * x]
    noiseeqs = [β * x]

    de = SDESystem(eqs, noiseeqs, t, sts, ps; name=:de)
    de = complete(de)

    x0 = 0.1
    u0map = [x => x0]

    parammap = [α => 1.5, β => 1.0]

    prob = SDEProblem(de, u0map, (0.0, 1.0), parammap)
    integrator = init(prob, LambaEulerHeun())
    # @show integrator |> typeof |> fieldnames
    solve(prob, LambaEulerHeun(); seed=1)

    # @show prob.p

    # prob.f.f.f_iip([1.0], [1.0], prob.p.tunable[1], 0.0)
end
