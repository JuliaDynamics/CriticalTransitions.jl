using CriticalTransitions
using BenchmarkTools

function benchmark_rate_system!(SUITE)
    function f(u, p, t) # out-of-place
        x = u[1]
        λ = p[1]
        dx = (x + λ)^2 - 1
        return SVector{1}(dx)
    end
    x0 = SVector{1}([-1.0])
    p_auto = [0.0]
    ds = CoupledODEs(f, x0, p_auto)

    p(t) = tanh(t)
    interval = (-100, 100)
    rc = RateFunction(p, interval)

    pidx = 1
    forcing_start_time = -100.0
    forcing_length = 200.0
    forcing_scale = 1.0
    t0 = forcing_start_time

    rs = RateSystem(ds, rc, pidx; forcing_start_time, forcing_length, forcing_scale, t0)

    T = forcing_length + 40.0
    trajectory(rs, T, x0)

    @btime trajectory($rs, $T, $x0)

    SUITE["Rate System"]["trajectory"]["RateSystem"] = @benchmarkable trajectory(
        $rs, $T, $x0
    ) seconds = 20

    function fexpl(u, p, t) # out-of-place
        x = u[1]
        λ = p[1]
        dx = (x + (tanh(t) + 1) / 2)^2 - 1
        return SVector{1}(dx)
    end
    x0 = SVector{1}([-1.0])
    p_auto = [0.0]
    t0 = -100
    nonauto_sysexpl = CoupledODEs(fexpl, x0, p_auto; t0)

    trajectory(nonauto_sysexpl, T, x0)
    @btime trajectory($nonauto_sysexpl, $T, $x0)
    return SUITE["Rate System"]["trajectory"]["Hard coded"] = @benchmarkable trajectory(
        $nonauto_sysexpl, $T, $x0
    ) seconds = 20
end
