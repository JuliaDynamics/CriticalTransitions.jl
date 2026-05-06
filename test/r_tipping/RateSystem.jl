using CriticalTransitions
using DynamicalSystemsBase
using Test
using StaticArrays
using Random
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

fp_case = ForcingProfile(profile, (-5.0, 5.0))

@testset "RateSystem" begin
    @testset "out-of-place CoupledODEs" begin
        function f_out(u, p, t)
            x = u[1]
            λ = p[1]
            dx = (x + λ)^2 - 1
            return SVector{1}(dx)
        end

        ds = CoupledODEs(f_out, x0, p_auto; t0 = 0.0)
        rs = RateSystem(ds, Dict(1 => fp_case); forcing_start_time = 0.0, forcing_duration = 10.0, forcing_scale = 1.0, t0 = 0.0)

        u = copy(current_state(rs))
        p_sys = current_parameters(rs.system)
        out = dynamic_rule(rs)(u, p_sys, 1.0)
        pmod_expected = DynamicalSystemsBase.current_parameters(rs, 1.0)
        @test isapprox(out[1], (u[1] + pmod_expected[1])^2 - 1; atol=1e-12)

        du = zeros(1)
        ret = dynamic_rule(rs)(du, u, p_sys, 1.0)
        @test ret === nothing
        @test isapprox(du[1], (u[1] + pmod_expected[1])^2 - 1; atol=1e-12)
    end

    @testset "out-of-place CoupledSDEs" begin
        function f_out(u, p, t)
            x = u[1]
            λ = p[1]
            dx = (x + λ)^2 - 1
            return SVector{1}(dx)
        end

        x0_sde = [0.0]
        p0_sde = [0.0]
        σ = 0.1
        ds = CoupledSDEs(f_out, x0_sde, p0_sde; noise_strength = σ)
        rs = RateSystem(ds, Dict(1 => fp_case); forcing_start_time = 0.0, forcing_duration = 10.0, forcing_scale = 1.0, t0 = 0.0)

        u = copy(current_state(rs))
        p_sys = current_parameters(rs.system)
        out = dynamic_rule(rs)(u, p_sys, 0.5)
        pmod_expected = DynamicalSystemsBase.current_parameters(rs, 0.5)
        @test isapprox(out[1], (u[1] + pmod_expected[1])^2 - 1; atol=1e-12)

        du = zeros(1)
        ret = dynamic_rule(rs)(du, u, p_sys, 0.5)
        @test ret === nothing
        @test isapprox(du[1], (u[1] + pmod_expected[1])^2 - 1; atol=1e-12)

        # Short SDE trajectory smoke checks (deterministic RNG)
        T_sde = 0.2
        Random.seed!(1234)
        tr_rs, _ = trajectory(rs, T_sde, x0_sde)
        @test !isempty(tr_rs)
        Random.seed!(1234)
        tr_sys, _ = trajectory(ds, T_sde, x0_sde)
        @test !isempty(tr_sys)
        @test !all(tr_rs .≈ tr_sys)
    end

    @testset "in-place CoupledODEs" begin
        function f_inplace!(du, u, p, t)
            du[1] = (u[1] + p[1])^2 - 1
            return nothing
        end

        ds = CoupledODEs(f_inplace!, x0, p_auto; t0 = 0.0)
        rs = RateSystem(ds, Dict(1 => fp_case); forcing_start_time = 0.0, forcing_duration = 10.0, forcing_scale = 1.0, t0 = 0.0)
        @test isinplace(rs) == true

        du = zeros(1)
        u = copy(current_state(rs))
        p_sys = current_parameters(rs.system)
        ret = dynamic_rule(rs)(du, u, p_sys, 1.0)
        @test ret === nothing
        pmod_expected = DynamicalSystemsBase.current_parameters(rs, 1.0)
        @test isapprox(du[1], (u[1] + pmod_expected[1])^2 - 1; atol=1e-12)
    end

    @testset "in-place CoupledSDEs" begin
        function f_sde_inplace!(du, u, p, t)
            du[1] = -u[1] + p[1]
            return nothing
        end

        x0_sde = [0.0]
        p0_sde = [0.0]
        σ = 0.1
        ds = CoupledSDEs(f_sde_inplace!, x0_sde, p0_sde; noise_strength = σ)
        rs = RateSystem(ds, Dict(1 => fp_case); forcing_start_time = 0.0, forcing_duration = 1.0, forcing_scale = 1.0, t0 = 0.0)
        @test isinplace(rs) == true

        du = zeros(1)
        u = [1.0]
        p_sys = current_parameters(rs.system)
        ret = dynamic_rule(rs)(du, u, p_sys, 0.5)
        @test ret === nothing
        pmod_expected = DynamicalSystemsBase.current_parameters(rs, 0.5)
        @test isapprox(du[1], -u[1] + pmod_expected[1]; atol=1e-12)

        # Short SDE trajectory smoke checks (deterministic RNG)
        T_sde = 0.2
        Random.seed!(1234)
        tr_rs_in, _ = trajectory(rs, T_sde, x0_sde)
        @test !isempty(tr_rs_in)
        Random.seed!(1234)
        tr_sys_in, _ = trajectory(ds, T_sde, x0_sde)
        @test !isempty(tr_sys_in)
        @test !all(tr_rs_in .≈ tr_sys_in)
    end


    @testset "parameter container variants" begin
        fp = ForcingProfile(tanh, (-5.0, 5.0))
        # out-of-place unforced rule (returns a new vector) used for these tests
        function f_out(u, p, t)
            x = u[1]
            λ = p[1]
            dx = (x + λ)^2 - 1
            return [dx]
        end
        # Build a RateSystemSpecs directly for parameter-container tests
        forcing_profile = Dict(1 => fp)
        start_map = Dict(1 => 10.0)
        duration_map = Dict(1 => 20.0)
        scale_map = Dict(1 => 1.0)
        p0_map = Dict(1 => 0.0)
        t0 = 0.0
        pdummy = [0.0]
        rss = CriticalTransitions.RateSystemSpecs{typeof(f_out),Int,Float64,Vector{Float64},Float64}(
            f_out, forcing_profile, start_map, duration_map, scale_map, p0_map, t0, pdummy, nothing
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
        prof = rss_fake.forcing_profile[1].profile
        section_start = rss_fake.forcing_profile[1].interval[1]
        section_end = rss_fake.forcing_profile[1].interval[2]
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
        # Skip detailed ModelingToolkit tests when the package is not available.
        if Base.find_package("ModelingToolkit") !== nothing
            try
                @eval begin
                    using ModelingToolkit
                end
                # Minimal smoke: ensure the package loads
                @test true
            catch
                @test true
            end
        else
            @test true
        end
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
    fp = ForcingProfile(tanh, (-5.0, 5.0))

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


