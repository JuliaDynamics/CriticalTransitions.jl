function _transition_step_loop(
        sys,
        x_i,
        x_f;
        radii = (0.1, 0.1),
        tmax = 1.0e3,
        radius_directions = 1:length(current_state(sys)),
        seed = nothing,
        kwargs...,
    )
    _, rad_f = radii
    prob = referenced_sciml_prob(sys)
    prob = CriticalTransitions.SciMLBase.remake(
        prob; u0 = oftype(prob.u0, x_i), tspan = (0, tmax)
    )
    if seed !== nothing
        prob = CriticalTransitions.SciMLBase.remake(prob; seed = UInt64(seed))
    end

    diffeq_kw = NamedTuple{filter(!=(:alg), keys(sys.diffeq))}(sys.diffeq)
    integ = CriticalTransitions.SciMLBase.init(
        prob, CriticalTransitions.solver(sys); diffeq_kw..., kwargs...
    )

    success = false
    while integ.t < tmax
        CriticalTransitions.SciMLBase.step!(integ)
        if CriticalTransitions.subnorm(integ.u - x_f; directions = radius_directions) < rad_f
            success = true
            CriticalTransitions.SciMLBase.terminate!(integ)
            break
        end
        CriticalTransitions.DynamicalSystemsBase.successful_step(integ) || break
    end

    sim = CriticalTransitions.SciMLBase.get_sol(integ)
    return StateSpaceSet(sim.u), sim.t, success
end

function benchmark_transition_callbacks!(suite)
    f(u, p, t) = [1.0]
    sys = CoupledSDEs(
        f,
        [0.0];
        noise_strength = 0.1,
        seed = 0x1234,
        diffeq = (
            alg = CriticalTransitions.StochasticDiffEq.EM(),
            dt = 1.0e-3,
            adaptive = false,
        ),
    )
    x_i = [0.0]
    x_f = [1.0]
    seed = 0x1234
    common = (; radii = (0.05, 0.05), tmax = 2.0, cut_start = false, seed)

    callback_result = transition(sys, x_i, x_f; common...)
    step_result = _transition_step_loop(
        sys,
        x_i,
        x_f;
        radii = common.radii,
        tmax = common.tmax,
        seed = common.seed,
    )
    @assert callback_result[3] == step_result[3]
    @assert callback_result[2] == step_result[2]
    @assert callback_result[1] == step_result[1]

    group = suite["Transitions"] = BenchmarkGroup()
    group["DiscreteCallback"] = @benchmarkable transition(
        $sys, $x_i, $x_f;
        radii = $(common.radii),
        tmax = $(common.tmax),
        cut_start = false,
        seed = $seed,
    )
    group["manual step loop"] = @benchmarkable _transition_step_loop(
        $sys, $x_i, $x_f;
        radii = $(common.radii),
        tmax = $(common.tmax),
        seed = $seed,
    )
    return suite
end
