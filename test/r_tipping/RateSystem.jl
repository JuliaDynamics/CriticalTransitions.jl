using CriticalTransitions
using DynamicalSystemsBase
using Test

using Random

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

mutable struct CustomParams
    λ::Float64
end
f_struct(u, p, t) = SVector((u[1] + p.λ)^2 - 1)

@testset "RateSystem" begin
    @testset "out-of-place CoupledODEs" begin
        ds = CoupledODEs(f, x0, p_auto; t0 = 0.0)
        rs = RateSystem(ds, Dict(1 => fp_case); forcing_start_time = 0.0, forcing_duration = 10.0, forcing_scale = 1.0, t0 = 0.0)

        u = copy(current_state(rs))
        p_sys = current_parameters(rs.system)
        out = dynamic_rule(rs)(u, p_sys, 1.0)
        pmod_expected = DynamicalSystemsBase.current_parameters(rs, 1.0)
        @test isapprox(out[1], (u[1] + pmod_expected[1])^2 - 1; atol = 1.0e-12)

        du = zeros(1)
        ret = dynamic_rule(rs)(du, u, p_sys, 1.0)
        @test ret === nothing
        @test isapprox(du[1], (u[1] + pmod_expected[1])^2 - 1; atol = 1.0e-12)
    end

    @testset "out-of-place CoupledSDEs" begin
        x0_sde = [0.0]
        p0_sde = [0.0]
        σ = 0.1
        ds = CoupledSDEs(f, x0_sde, p0_sde; noise_strength = σ)
        rs = RateSystem(ds, Dict(1 => fp_case); forcing_start_time = 0.0, forcing_duration = 10.0, forcing_scale = 1.0, t0 = 0.0)

        u = copy(current_state(rs))
        p_sys = current_parameters(rs.system)
        out = dynamic_rule(rs)(u, p_sys, 0.5)
        pmod_expected = DynamicalSystemsBase.current_parameters(rs, 0.5)
        @test isapprox(out[1], (u[1] + pmod_expected[1])^2 - 1; atol = 1.0e-12)

        du = zeros(1)
        ret = dynamic_rule(rs)(du, u, p_sys, 0.5)
        @test ret === nothing
        @test isapprox(du[1], (u[1] + pmod_expected[1])^2 - 1; atol = 1.0e-12)

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
        @test isapprox(du[1], (u[1] + pmod_expected[1])^2 - 1; atol = 1.0e-12)
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
        @test isapprox(du[1], -u[1] + pmod_expected[1]; atol = 1.0e-12)

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

    @testset "out-of-place CoupledODEs_2param" begin
        function f_out(u, p, t)
            x = u[1]
            λ = p[1]
            γ = p[2]
            dx = (x + λ + γ)^2 - 1
            return SVector{1}(dx)
        end

        ds = CoupledODEs(f_out, x0, [0.0, 0.0]; t0 = 0.0)
        rs = RateSystem(
            ds,
            Dict(1 => fp_case, 2 => fp_case);
            forcing_start_time = Dict(1 => 0.0, 2 => 0.0),
            forcing_duration = Dict(1 => 10.0, 2 => 10.0),
            forcing_scale = Dict(1 => 0.5, 2 => 0.5),
            t0 = 0.0,
        )

        u = copy(current_state(rs))
        p_sys = current_parameters(rs.system)
        out = dynamic_rule(rs)(u, p_sys, 1.0)
        pmod_expected = DynamicalSystemsBase.current_parameters(rs, 1.0)
        peff = pmod_expected[1] + pmod_expected[2]
        @test isapprox(out[1], (u[1] + peff)^2 - 1; atol = 1.0e-12)

        du = zeros(1)
        ret = dynamic_rule(rs)(du, u, p_sys, 1.0)
        @test ret === nothing
        @test isapprox(du[1], (u[1] + peff)^2 - 1; atol = 1.0e-12)
    end

    @testset "out-of-place CoupledSDEs_2param" begin
        function f_out(u, p, t)
            x = u[1]
            λ = p[1]
            γ = p[2]
            dx = (x + λ + γ)^2 - 1
            return SVector{1}(dx)
        end

        x0_sde = [0.0]
        p0_sde = [0.0, 0.0]
        σ = 0.1
        ds = CoupledSDEs(f_out, x0_sde, p0_sde; noise_strength = σ)
        rs = RateSystem(
            ds,
            Dict(1 => fp_case, 2 => fp_case);
            forcing_start_time = Dict(1 => 0.0, 2 => 0.0),
            forcing_duration = Dict(1 => 10.0, 2 => 10.0),
            forcing_scale = Dict(1 => 0.5, 2 => 0.5),
            t0 = 0.0,
        )

        u = copy(current_state(rs))
        p_sys = current_parameters(rs.system)
        out = dynamic_rule(rs)(u, p_sys, 0.5)
        pmod_expected = DynamicalSystemsBase.current_parameters(rs, 0.5)
        peff = pmod_expected[1] + pmod_expected[2]
        @test isapprox(out[1], (u[1] + peff)^2 - 1; atol = 1.0e-12)

        du = zeros(1)
        ret = dynamic_rule(rs)(du, u, p_sys, 0.5)
        @test ret === nothing
        @test isapprox(du[1], (u[1] + peff)^2 - 1; atol = 1.0e-12)

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

    @testset "reverse forcing" begin
        ds = CoupledODEs(f, x0, p_auto; t0 = 0.0)
        rs = RateSystem(
            ds,
            Dict(1 => fp_case);
            forcing_start_time = 0.0,
            forcing_duration = 10.0,
            forcing_scale = 1.0,
            reverse = true,
            t0 = 0.0,
        )

        f0 = fp_case.profile(fp_case.interval[1])
        fend = fp_case.profile(fp_case.interval[2])

        # End of forward interval: reaches the same value as non-reversed forcing.
        @test isapprox(parameter(rs, 10.0, 1), fend - f0; atol = 1.0e-12)

        # During second interval: follows profile backwards.
        reverse_mid_t = 15.0
        reverse_shift = fp_case.interval[2] - ((fp_case.interval[2] - fp_case.interval[1]) / 10.0) * (reverse_mid_t - 10.0)
        expected_mid = fp_case.profile(reverse_shift) - f0
        @test isapprox(parameter(rs, reverse_mid_t, 1), expected_mid; atol = 1.0e-12)

        # After reverse phase ends: returns to autonomous value.
        @test isapprox(parameter(rs, 21.0, 1), 0.0; atol = 1.0e-12)
    end

    @testset "reverse as dictionary" begin
        function f_out(u, p, t)
            x = u[1]
            λ = p[1]
            γ = p[2]
            dx = (x + λ + γ)^2 - 1
            return SVector{1}(dx)
        end

        ds = CoupledODEs(f_out, x0, [0.0, 0.0]; t0 = 0.0)
        rs = RateSystem(
            ds,
            Dict(1 => fp_case, 2 => fp_case);
            forcing_start_time = Dict(1 => 0.0, 2 => 0.0),
            forcing_duration = Dict(1 => 10.0, 2 => 10.0),
            forcing_scale = Dict(1 => 1.0, 2 => 1.0),
            reverse = Dict(1 => true, 2 => false),
            t0 = 0.0,
        )

        f0 = fp_case.profile(fp_case.interval[1])
        fend = fp_case.profile(fp_case.interval[2])

        # First parameter reversed: returns to base after 2*duration.
        @test isapprox(parameter(rs, 21.0, 1), 0.0; atol = 1.0e-12)
        # Second parameter non-reversed: stays at final forced value.
        @test isapprox(parameter(rs, 21.0, 2), fend - f0; atol = 1.0e-12)
    end


    @testset "custom (non-Vector) parameter container" begin
        # Forcing should work for any parameter container addressable by `set_parameter!`,
        # here a mutable struct keyed by the Symbol `:λ`, not just `Vector`s.
        ds = CoupledODEs(f_struct, x0, CustomParams(0.0); t0 = 0.0)
        start_time, duration, scale = 10.0, 20.0, 1.0
        rs = RateSystem(
            ds, fp_case, :λ;
            forcing_start_time = start_time, forcing_duration = duration,
            forcing_scale = scale, t0 = 0.0,
        )

        section_start, section_end = fp_case.interval
        # Mid-forcing: parameter follows the rescaled profile.
        t = 15.0
        time_shift = ((section_end - section_start) / duration) * (t - start_time) + section_start
        expected = scale * (fp_case.profile(time_shift) - fp_case.profile(section_start))
        @test isapprox(parameter(rs, t, :λ), expected; atol = 1.0e-12)

        # Before forcing starts: parameter at its autonomous value.
        @test isapprox(parameter(rs, 5.0, :λ), 0.0; atol = 1.0e-12)
        # After forcing ends: frozen at the final value.
        fend = scale * (fp_case.profile(section_end) - fp_case.profile(section_start))
        @test isapprox(parameter(rs, start_time + duration + 1.0, :λ), fend; atol = 1.0e-12)

        # The forcing must drive the actual rule and preserve the container type.
        @test current_parameters(rs.system) isa CustomParams
        u = [2.0]
        out = dynamic_rule(rs)(u, current_parameters(rs.system), t)
        @test isapprox(out[1], (u[1] + parameter(rs, t, :λ))^2 - 1; atol = 1.0e-12)
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

    rs = RateSystem(ds, Dict(pidx => fp); forcing_start_time, forcing_duration, forcing_scale)

    @testset "unforced_system" begin
        frozen = unforced_system(rs, rs.specs.t0)
        unforced_rule = dynamic_rule(rs.specs.unforced_system)
        u = [2.0]
        t = 10.0
        du_frozen = dynamic_rule(frozen)(u, current_parameters(frozen), t)
        du_orig = dynamic_rule(ds)(u, current_parameters(ds), t)
        du_unforced = unforced_rule(u, initial_parameters(rs), t)

        @test du_frozen == du_orig
        @test du_frozen == du_unforced
    end

    @testset "unforced_system preserves CoupledSDEs" begin
        ds_sde = CoupledSDEs(f, [0.0], [0.0]; noise_strength = 0.1)
        rs_sde = RateSystem(
            ds_sde, Dict(1 => fp_case);
            forcing_start_time = 0.0, forcing_duration = 10.0, forcing_scale = 1.0, t0 = 0.0
        )
        unforced_sde = unforced_system(rs_sde, 0.0)
        @test unforced_sde isa CoupledSDEs
        @test !isdeterministic(unforced_sde)
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
            ds, Dict(pidx => fp); forcing_start_time, forcing_duration, forcing_scale
        )

        p_hardcoded = [p_auto[1], forcing_scale, forcing_duration / 2, 40 / forcing_duration]
        sys_hardcoded = CoupledODEs(fexpl, x0, p_hardcoded; t0 = 0.0)

        tr_constructed, _ = trajectory(sys_constructed, T, x0)
        tr_hardcoded, _ = trajectory(sys_hardcoded, T, x0)

        @test all(tr_constructed .≈ tr_hardcoded)
    end
end
