function _linear_transition_case(dim, radius_directions)
    f!(du, u, p, t) = (fill!(du, 1.0); nothing)
    sys = CoupledSDEs(
        f!,
        zeros(dim);
        noise_strength = 0.02,
        seed = 0x1234,
        diffeq = (
            alg = CriticalTransitions.StochasticDiffEq.EM(),
            dt = 1.0e-3,
            adaptive = false,
        ),
    )
    x_i = zeros(dim)
    x_f = ones(dim)
    rad = 0.1 * sqrt(length(radius_directions))
    seed = 0x1234
    common = (; radii = (rad, rad), tmax = 2.0, cut_start = false, seed, radius_directions)

    original = transition(sys, x_i, x_f; common...)
    optimized = _transition_callback_noalloc(
        sys,
        x_i,
        x_f;
        radii = common.radii,
        tmax = common.tmax,
        radius_directions = common.radius_directions,
        seed = common.seed,
    )
    @assert original[3] == optimized[3]
    @assert original[2] == optimized[2]
    @assert original[1] == optimized[1]

    return (; sys, x_i, x_f, common)
end

function _add_transition_condition_case!(group, name, dim, radius_directions)
    case = _linear_transition_case(dim, radius_directions)
    sys = case.sys
    x_i = case.x_i
    x_f = case.x_f
    common = case.common

    subgroup = group[name] = BenchmarkGroup()
    subgroup["original subnorm"] = @benchmarkable transition(
        $sys,
        $x_i,
        $x_f;
        radii = $(common.radii),
        tmax = $(common.tmax),
        cut_start = false,
        radius_directions = $(common.radius_directions),
        seed = $(common.seed),
    )
    subgroup["squared distance"] = @benchmarkable _transition_callback_noalloc(
        $sys,
        $x_i,
        $x_f;
        radii = $(common.radii),
        tmax = $(common.tmax),
        radius_directions = $(common.radius_directions),
        seed = $(common.seed),
    )
    return group
end

function benchmark_transition_condition_scaling!(suite)
    group = suite["Transition target condition scaling"] = BenchmarkGroup()
    _add_transition_condition_case!(group, "D=1, all", 1, 1:1)
    _add_transition_condition_case!(group, "D=2, all", 2, 1:2)
    _add_transition_condition_case!(group, "D=8, all", 8, 1:8)
    _add_transition_condition_case!(group, "D=8, first 2", 8, 1:2)
    _add_transition_condition_case!(group, "D=32, all", 32, 1:32)
    _add_transition_condition_case!(group, "D=32, sparse 5", 32, [1, 8, 16, 24, 32])
    return suite
end
