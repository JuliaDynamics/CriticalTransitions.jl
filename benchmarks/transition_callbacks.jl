function _inside_target_noalloc(u, x_f, rad_f, radius_directions)
    d2 = zero(promote_type(eltype(u), eltype(x_f), typeof(rad_f)))
    @inbounds for i in radius_directions
        δ = u[i] - x_f[i]
        d2 += δ * δ
    end
    return d2 < rad_f * rad_f
end

function _transition_step_loop(
        sys,
        x_i,
        x_f;
        radii = (0.1, 0.1),
        tmax = 1.0e3,
        radius_directions = 1:length(current_state(sys)),
        seed = nothing,
        check_success = true,
        allocation_free_condition = false,
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
        reached = if allocation_free_condition
            _inside_target_noalloc(integ.u, x_f, rad_f, radius_directions)
        else
            CriticalTransitions.subnorm(integ.u - x_f; directions = radius_directions) < rad_f
        end
        if reached
            success = true
            CriticalTransitions.SciMLBase.terminate!(integ)
            CriticalTransitions.SciMLBase.savevalues!(integ, true)
            break
        end
        check_success && !CriticalTransitions.DynamicalSystemsBase.successful_step(integ) && break
    end

    sim = CriticalTransitions.SciMLBase.get_sol(integ)
    return StateSpaceSet(sim.u), sim.t, success
end

function _run_transition_steps!(integ, x_f, rad_f, radius_directions, tmax)
    success = false
    while integ.t < tmax
        CriticalTransitions.SciMLBase.step!(integ)
        if _inside_target_noalloc(integ.u, x_f, rad_f, radius_directions)
            success = true
            CriticalTransitions.SciMLBase.terminate!(integ)
            CriticalTransitions.SciMLBase.savevalues!(integ, true)
            break
        end
        CriticalTransitions.DynamicalSystemsBase.successful_step(integ) || break
    end
    return success
end

function _transition_step_loop_barrier(
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
    success = _run_transition_steps!(integ, x_f, rad_f, radius_directions, tmax)
    sim = CriticalTransitions.SciMLBase.get_sol(integ)
    return StateSpaceSet(sim.u), sim.t, success
end

function _transition_callback_noalloc(
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
    condition(u, t, integrator) = _inside_target_noalloc(u, x_f, rad_f, radius_directions)
    affect!(integrator) = CriticalTransitions.SciMLBase.terminate!(integrator)
    cb_ball = CriticalTransitions.SciMLBase.DiscreteCallback(condition, affect!)

    prob = referenced_sciml_prob(sys)
    prob = CriticalTransitions.SciMLBase.remake(
        prob; u0 = oftype(prob.u0, x_i), tspan = (0, tmax)
    )
    if seed !== nothing
        prob = CriticalTransitions.SciMLBase.remake(prob; seed = UInt64(seed))
    end
    diffeq_kw = NamedTuple{filter(!=(:alg), keys(sys.diffeq))}(sys.diffeq)
    sim = CriticalTransitions.SciMLBase.solve(
        prob, CriticalTransitions.solver(sys); callback = cb_ball, diffeq_kw..., kwargs...
    )
    success = sim.retcode == CriticalTransitions.SciMLBase.ReturnCode.Terminated
    return StateSpaceSet(sim.u), sim.t, success
end

function _solve_horizon(sys, x_i, tend; seed = nothing, kwargs...)
    prob = referenced_sciml_prob(sys)
    prob = CriticalTransitions.SciMLBase.remake(
        prob; u0 = oftype(prob.u0, x_i), tspan = (0, tend)
    )
    if seed !== nothing
        prob = CriticalTransitions.SciMLBase.remake(prob; seed = UInt64(seed))
    end
    diffeq_kw = NamedTuple{filter(!=(:alg), keys(sys.diffeq))}(sys.diffeq)
    return CriticalTransitions.SciMLBase.solve(
        prob, CriticalTransitions.solver(sys); diffeq_kw..., kwargs...
    )
end

function _step_horizon(sys, x_i, tend; seed = nothing, kwargs...)
    prob = referenced_sciml_prob(sys)
    prob = CriticalTransitions.SciMLBase.remake(
        prob; u0 = oftype(prob.u0, x_i), tspan = (0, tend)
    )
    if seed !== nothing
        prob = CriticalTransitions.SciMLBase.remake(prob; seed = UInt64(seed))
    end
    diffeq_kw = NamedTuple{filter(!=(:alg), keys(sys.diffeq))}(sys.diffeq)
    integ = CriticalTransitions.SciMLBase.init(
        prob, CriticalTransitions.solver(sys); diffeq_kw..., kwargs...
    )
    while integ.t < tend
        CriticalTransitions.SciMLBase.step!(integ)
    end
    return CriticalTransitions.SciMLBase.get_sol(integ)
end

function _run_horizon_steps!(integ, tend)
    while integ.t < tend
        CriticalTransitions.SciMLBase.step!(integ)
    end
    return integ
end

function _step_horizon_barrier(sys, x_i, tend; seed = nothing, kwargs...)
    prob = referenced_sciml_prob(sys)
    prob = CriticalTransitions.SciMLBase.remake(
        prob; u0 = oftype(prob.u0, x_i), tspan = (0, tend)
    )
    if seed !== nothing
        prob = CriticalTransitions.SciMLBase.remake(prob; seed = UInt64(seed))
    end
    diffeq_kw = NamedTuple{filter(!=(:alg), keys(sys.diffeq))}(sys.diffeq)
    integ = CriticalTransitions.SciMLBase.init(
        prob, CriticalTransitions.solver(sys); diffeq_kw..., kwargs...
    )
    _run_horizon_steps!(integ, tend)
    return CriticalTransitions.SciMLBase.get_sol(integ)
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
    noalloc_callback_result = _transition_callback_noalloc(
        sys,
        x_i,
        x_f;
        radii = common.radii,
        tmax = common.tmax,
        seed = common.seed,
    )
    noalloc_step_result = _transition_step_loop(
        sys,
        x_i,
        x_f;
        radii = common.radii,
        tmax = common.tmax,
        seed = common.seed,
        allocation_free_condition = true,
    )
    barrier_step_result = _transition_step_loop_barrier(
        sys,
        x_i,
        x_f;
        radii = common.radii,
        tmax = common.tmax,
        seed = common.seed,
    )
    @assert callback_result[3] == step_result[3] == noalloc_callback_result[3] == noalloc_step_result[3] == barrier_step_result[3]
    @assert callback_result[2] == step_result[2] == noalloc_callback_result[2] == noalloc_step_result[2] == barrier_step_result[2]
    @assert callback_result[1] == step_result[1] == noalloc_callback_result[1] == noalloc_step_result[1] == barrier_step_result[1]

    terminal_time = callback_result[2][end]
    solve_horizon = _solve_horizon(sys, x_i, terminal_time; seed)
    step_horizon = _step_horizon(sys, x_i, terminal_time; seed)
    barrier_horizon = _step_horizon_barrier(sys, x_i, terminal_time; seed)
    @assert solve_horizon.t == step_horizon.t == barrier_horizon.t
    @assert solve_horizon.u == step_horizon.u == barrier_horizon.u
    @info "transition benchmark" saved_points = length(callback_result[2]) terminal_time

    group = suite["Transitions"] = BenchmarkGroup()
    group["DiscreteCallback"] = @benchmarkable transition(
        $sys, $x_i, $x_f;
        radii = $(common.radii),
        tmax = $(common.tmax),
        cut_start = false,
        seed = $seed,
    )
    group["DiscreteCallback noalloc condition"] = @benchmarkable _transition_callback_noalloc(
        $sys, $x_i, $x_f;
        radii = $(common.radii),
        tmax = $(common.tmax),
        seed = $seed,
    )
    group["manual step loop"] = @benchmarkable _transition_step_loop(
        $sys, $x_i, $x_f;
        radii = $(common.radii),
        tmax = $(common.tmax),
        seed = $seed,
    )
    group["manual step loop no success check"] = @benchmarkable _transition_step_loop(
        $sys, $x_i, $x_f;
        radii = $(common.radii),
        tmax = $(common.tmax),
        seed = $seed,
        check_success = false,
    )
    group["manual step loop noalloc condition"] = @benchmarkable _transition_step_loop(
        $sys, $x_i, $x_f;
        radii = $(common.radii),
        tmax = $(common.tmax),
        seed = $seed,
        allocation_free_condition = true,
    )
    group["manual step loop noalloc barrier"] = @benchmarkable _transition_step_loop_barrier(
        $sys, $x_i, $x_f;
        radii = $(common.radii),
        tmax = $(common.tmax),
        seed = $seed,
    )
    group["solve fixed horizon no callback"] = @benchmarkable _solve_horizon(
        $sys, $x_i, $terminal_time; seed = $seed
    )
    group["step fixed horizon no callback"] = @benchmarkable _step_horizon(
        $sys, $x_i, $terminal_time; seed = $seed
    )
    group["step fixed horizon barrier"] = @benchmarkable _step_horizon_barrier(
        $sys, $x_i, $terminal_time; seed = $seed
    )
    return suite
end
