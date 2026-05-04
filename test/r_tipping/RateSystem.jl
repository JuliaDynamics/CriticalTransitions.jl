using CriticalTransitions
using DynamicalSystemsBase
using Test
using StaticArrays
using ModelingToolkit

pidx = 1
forcing_start_time = 20.0
forcing_duration = 100.0
forcing_scale = 1.5
T = forcing_duration + 40.0
x0 = [-1.0]
p_auto = [0.0]

profile(t) = tanh(t)

# Autonomous drift
function f(u, p, t) # out-of-place
    x = u[1]
    λ = p[1]
    dx = (x + λ)^2 - 1
    return SVector{1}(dx)
end


@testset "call-signature compatibility" begin
    # In-place unforced rule (mutates `du`) ---------------------------------
    function f_inplace!(du, u, p, t)
        du[1] = (u[1] + p[1])^2 - 1
        return nothing
    end

    # Construct RateSystemSpecs directly to avoid constructing CoupledODEs
    fp_ip = ForcingProfile(profile, (-5.0, 5.0))
    forcers_ip = Dict(1 => fp_ip)
    start_map_ip = Dict(1 => 10.0)
    duration_map_ip = Dict(1 => 20.0)
    scale_map_ip = Dict(1 => 1.0)
    p0_map_ip = Dict(1 => 0.0)
    t0_ip = 0.0
    pdummy_ip = [0.0]
    rss_ip = CriticalTransitions.RateSystemSpecs{typeof(f_inplace!),Int,Float64,Vector{Float64},Float64}(
        f_inplace!, forcers_ip, start_map_ip, duration_map_ip, scale_map_ip, p0_map_ip, t0_ip, pdummy_ip, nothing
    )

    u = [2.0]
    du = zeros(1)
    pcopy = deepcopy(rss_ip.pdummy)
    t_ip = 5.0
    pmod_expected = CriticalTransitions.p_modified(rss_ip, deepcopy(pcopy), t_ip)
    expected = (u[1] + pmod_expected[1])^2 - 1

    ret = rss_ip(du, u, pcopy, t_ip)
    @test ret === nothing
    @test isapprox(du[1], expected; atol=1e-12)

    # Out-of-place unforced rule (returns a new vector) ---------------------
    function f_out(u, p, t)
        x = u[1]
        λ = p[1]
        dx = (x + λ)^2 - 1
        return [dx]
    end

    # out-of-place: RateSystemSpecs with Vector parameter container
    fp_op = ForcingProfile(profile, (-5.0, 5.0))
    forcers_op = Dict(1 => fp_op)
    start_map_op = Dict(1 => 10.0)
    duration_map_op = Dict(1 => 20.0)
    scale_map_op = Dict(1 => 1.0)
    p0_map_op = Dict(1 => 0.0)
    t0_op = 0.0
    pdummy_op = [0.0]
    rss_op = CriticalTransitions.RateSystemSpecs{typeof(f_out),Int,Float64,Vector{Float64},Float64}(
        f_out, forcers_op, start_map_op, duration_map_op, scale_map_op, p0_map_op, t0_op, pdummy_op, nothing
    )

    u2 = [2.0]
    du2 = zeros(1)
    pcopy2 = deepcopy(rss_op.pdummy)
    t_op = 5.0
    pmod_expected2 = CriticalTransitions.p_modified(rss_op, deepcopy(pcopy2), t_op)
    expected2 = (u2[1] + pmod_expected2[1])^2 - 1

    ret2 = rss_op(du2, u2, pcopy2, t_op)
    @test ret2 === nothing
    @test isapprox(du2[1], expected2; atol=1e-12)
end


@testset "parameter container variants" begin
    fp = ForcingProfile(profile, (-5.0, 5.0))
    # out-of-place unforced rule (returns a new vector) used for these tests
    function f_out(u, p, t)
        x = u[1]
        λ = p[1]
        dx = (x + λ)^2 - 1
        return [dx]
    end
    # Build a RateSystemSpecs directly for parameter-container tests
    forcers = Dict(1 => fp)
    start_map = Dict(1 => 10.0)
    duration_map = Dict(1 => 20.0)
    scale_map = Dict(1 => 1.0)
    p0_map = Dict(1 => 0.0)
    t0 = 0.0
    pdummy = [0.0]
    rss = CriticalTransitions.RateSystemSpecs{typeof(f_out),Int,Float64,Vector{Float64},Float64}(
        f_out, forcers, start_map, duration_map, scale_map, p0_map, t0, pdummy, nothing
    )
    t = 5.0

    # Dict: should be mutated in-place and returned
    pd = Dict(1 => 0.0)
    pd_ret = CriticalTransitions.p_modified(rss, pd, t)
    @test pd_ret === pd

    # Vector: should be mutated in-place and returned
    pv = [0.0]
    pv_ret = CriticalTransitions.p_modified(rss, pv, t)
    @test pv_ret === pv

    # Arbitrary container via owner/set_parameter! --------------------------------
    mutable struct FakeOwner
        params::Dict{Int, Float64}
    end

    function DynamicalSystemsBase.set_parameter!(o::FakeOwner, k, v)
        o.params[k] = v
        return o
    end

    function DynamicalSystemsBase.current_parameters(o::FakeOwner)
        return o.params
    end

    function DynamicalSystemsBase.set_parameters!(o::FakeOwner, pd)
        o.params = deepcopy(pd)
        return o
    end

    # construct a RateSystemSpecs and then attach the fake owner/pdummy
    forcers_fake = Dict(1 => fp)
    start_map_fake = Dict(1 => 0.0)
    duration_map_fake = Dict(1 => 10.0)
    scale_map_fake = Dict(1 => 1.0)
    p0_map_fake = Dict(1 => 0.0)
    t0_fake = 0.0
    pdummy_fake = Dict(1 => 0.0)
    rss_fake = CriticalTransitions.RateSystemSpecs{typeof(f_out),Int,Float64,Any,Float64}(
        f_out, forcers_fake, start_map_fake, duration_map_fake, scale_map_fake, p0_map_fake, t0_fake, pdummy_fake, nothing
    )
    fake = FakeOwner(Dict(1 => 0.0))
    rss_fake.owner = fake
    rss_fake.pdummy = deepcopy(fake.params)

    pd_fake = CriticalTransitions.p_modified(rss_fake, Dict(1 => 0.0), 5.0)
    # compute expected updated parameter value using the configured profile
    p0 = initial_parameter(rss_fake, 1)
    prof = rss_fake.forcers[1].profile
    section_start = rss_fake.forcers[1].interval[1]
    section_end = rss_fake.forcers[1].interval[2]
    start_time = rss_fake.forcing_start_time[1]
    duration = rss_fake.forcing_duration[1]
    scale = rss_fake.forcing_scale[1]

    if 5.0 <= start_time
        expected_pt = p0
    elseif 5.0 < start_time + duration
        time_shift = ((section_end - section_start) / duration) * (5.0 - start_time) + section_start
        expected_pt = p0 + scale * (prof(time_shift) - prof(section_start))
    else
        expected_pt = p0 + scale * (prof(section_end) - prof(section_start))
    end

    @test haskey(pd_fake, 1)
    @test isapprox(pd_fake[1], expected_pt; atol=1e-12)
end


@testset "ModelingToolkit smoke" begin
    # Basic ModelingToolkit construction smoke-check (no solve)
    @independent_variables t
    @parameters a
    @variables x(t)
    D = Differential(t)
    eqs = [D(x) ~ (x + a)^2 - 1]
    @named sys = ODESystem(eqs, t)
    @test length(eqs) == 1
    @test sys !== nothing
end

# Hard-coded non-autonomous drift
function fexpl(u, p, t) # out-of-place
    x = u[1]
    λ_0, scale, shift, rate = p
    λ = λ_0 + scale * (tanh(rate * (t - shift)) + 1)
    dx = (x + λ)^2 - 1
    return SVector{1}(dx)
end

ds = CoupledODEs(f, x0, p_auto)

profile(t) = tanh(t)
fp = ForcingProfile(profile, (-5.0, 5.0))

rs = RateSystem(ds, Dict(pidx => fp); forcing_start_time = forcing_start_time, forcing_duration = forcing_duration, forcing_scale = forcing_scale)

@testset "frozen_system" begin
    frozen = frozen_system(rs, rs.forcing.t0)
    unforced = rs.forcing.unforced_rule
    u = [2.0]
    t = 10.0
    du_frozen = dynamic_rule(frozen)(u, current_parameters(frozen), t)
    du_orig = dynamic_rule(ds)(u, current_parameters(ds), t)
    du_unforced = unforced(u, [initial_parameter_value(rs)], t)

    @test du_frozen == du_orig
    @test du_frozen == du_unforced
end

@testset "RateSystem" begin
    @testset "DynamicalSystems API" begin
        @test current_state(rs) == x0
        @test DynamicalSystemsBase.initial_state(rs) == x0
        @test current_parameters(rs) == p_auto
        @test initial_parameters(rs) == p_auto
        @test current_time(rs) == 0.0
        @test initial_time(rs) == 0.0
        @test isinplace(rs) == false
        @test isdeterministic(rs) == true
        @test isdiscretetime(rs) == false
        @test rs(0) == x0
    end

    # Compute trajectories
    auto_traj = trajectory(ds, T, x0)
    nonauto_traj = trajectory(rs, T, x0)

    @test isapprox(auto_traj[1][end, 1], -1; atol = 1.0e-2)
    @test isapprox(nonauto_traj[1][end, 1], -4; atol = 1.0e-2)

    @testset "Equivalence with hard-coded" begin
        forcing_start_time = 0.0
        forcing_duration = 100.0
        forcing_scale = 1.0
        fp = ForcingProfile(profile, (-20.0, 20.0))

        ds = CoupledODEs(f, x0, p_auto)
        sys_constructed = RateSystem(
            ds, Dict(pidx => fp); forcing_start_time = forcing_start_time, forcing_duration = forcing_duration, forcing_scale = forcing_scale
        )

        p_hardcoded = [p_auto[1], forcing_scale, forcing_duration / 2, 40 / forcing_duration]
        sys_hardcoded = CoupledODEs(fexpl, x0, p_hardcoded; t0 = 0.0)

        tr_constructed, _ = trajectory(sys_constructed, T, x0)
        tr_hardcoded, _ = trajectory(sys_hardcoded, T, x0)

        @test all(tr_constructed .≈ tr_hardcoded)
    end
end

@testset "RateSystem SDE smoke" begin
    # simple scalar SDE where the parameter p[1] enters the drift
    function f_sde(u, p, t)
        return [ -u[1] + p[1] ]
    end

    x0_sde = [0.0]
    p0_sde = [0.0]
    σ = 0.1

    sys = CoupledSDEs(f_sde, x0_sde, p0_sde; noise_strength=σ)

    fp = ForcingProfile(:linear)
    rs_sde = RateSystem(sys, Dict(1 => fp); forcing_start_time=0.0, forcing_duration=0.1, forcing_scale=0.5, t0=0.0)

    T = 0.2
    tr_rs, _ = trajectory(rs_sde, T, x0_sde)
    @test !isempty(tr_rs)

    # Also compute a short trajectory of the underlying SDE directly
    tr_sys, _ = trajectory(sys, T, x0_sde)
    @test !isempty(tr_sys)
end
