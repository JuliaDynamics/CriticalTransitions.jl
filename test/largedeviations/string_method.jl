using CriticalTransitions
using OrdinaryDiffEq: Tsit5
using LinearAlgebra: norm

const λ = 3 / 1.21 * 2 / 295
const ω0 = 1.000
const ω = 1.000
const γ = 1 / 295
const η = 0
const α = -1

function fu(u, v)
    return (-4 * γ * ω * u - 2 * λ * v - 4 * (ω0 - ω^2) * v - 3 * α * v * (u^2 + v^2)) /
           (8 * ω)
end
function fv(u, v)
    return (-4 * γ * ω * v - 2 * λ * u + 4 * (ω0 - ω^2) * u + 3 * α * u * (u^2 + v^2)) /
           (8 * ω)
end
stream(u, v) = Point2f(fu(u, v), fv(u, v))
dfvdv(u, v) = (-4 * γ * ω + 6 * α * u * v) / (8 * ω)
dfudu(u, v) = (-4 * γ * ω - 6 * α * u * v) / (8 * ω)
dfvdu(u, v) = (-2 * λ + 4 * (ω0 - ω^2) + 9 * α * u^2 + 3 * α * v^2) / (8 * ω)
dfudv(u, v) = (-2 * λ - 4 * (ω0 - ω^2) - 3 * α * u^2 - 9 * α * v^2) / (8 * ω)

Nt = 500  # number of discrete time steps
s = collect(range(0; stop=1, length=Nt))

xa = [-0.0208, 0.0991]
xb = -xa
xsaddle = [0.0, 0.0]

# Initial trajectory
xx = @. (xb[1] - xa[1]) * s + xa[1] + 4 * s * (1 - s) * xsaddle[1]
yy = @. (xb[2] - xa[2]) * s + xa[2] + 4 * s * (1 - s) * xsaddle[2] + 0.01 * sin(2π * s)

@testset "StateSpaceSet vs Matrix" begin
    function H_x_sss(x, p) # ℜ² → ℜ²
        u, v = eachcol(x)
        pu, pv = eachcol(p)

        H_u = @. pu * dfudu(u, v) + pv * dfvdu(u, v)
        H_v = @. pu * dfudv(u, v) + pv * dfvdv(u, v)
        return StateSpaceSet([H_u H_v])
    end
    function H_p_sss(x, p) # ℜ² → ℜ²
        u, v = eachcol(x)
        pu, pv = eachcol(p)

        H_pu = @. pu + fu(u, v)
        H_pv = @. pv + fv(u, v)
        return StateSpaceSet([H_pu H_pv])
    end

    sys_sss = ExtendedPhaseSpace{false,2}(H_x_sss, H_p_sss)

    x_init_sss = StateSpaceSet([xx yy])

    string_sss = string_method(
        sys_sss, x_init_sss; iterations=10_000, ϵ=0.5, show_progress=false
    )

    function H_x_m(x, p) # ℜ² → ℜ²
        u, v = eachrow(x)
        pu, pv = eachrow(p)

        H_u = @. pu * dfudu(u, v) + pv * dfvdu(u, v)
        H_v = @. pu * dfudv(u, v) + pv * dfvdv(u, v)
        return Matrix([H_u H_v]')
    end
    function H_p_m(x, p) # ℜ² → ℜ²
        u, v = eachrow(x)
        pu, pv = eachrow(p)

        H_pu = @. pu + fu(u, v)
        H_pv = @. pv + fv(u, v)
        return Matrix([H_pu H_pv]')
    end

    sys_m = ExtendedPhaseSpace{false,2}(H_x_m, H_p_m)

    x_init_m = Matrix([xx yy]')

    string_m = string_method(sys_m, x_init_m; iterations=10_000, ϵ=0.5, show_progress=false)

    @test vec(string_m) ≈ vec(string_sss)
end

@testset "String method endpoints pinned" begin
    Nt_small = 50
    s_small = range(0; stop=1, length=Nt_small)

    xa2 = [-1.0, 0.0]
    xb2 = [1.0, 0.0]

    x_init_m =
        xa2 .* (1 .- s_small)' .+ xb2 .* s_small' .+ [0.0, 0.3] .* sinpi.(s_small)'

    b_rot(x) = [-x[2], x[1]]

    string_m = string_method(
        b_rot,
        x_init_m;
        iterations=25,
        ϵ=0.2,
        show_progress=false,
    )
    m = Matrix(string_m)
    @test vec(m[1, :]) ≈ x_init_m[:, 1]
    @test vec(m[end, :]) ≈ x_init_m[:, end]

    x_init_sss = StateSpaceSet(x_init_m')
    string_sss = string_method(
        b_rot,
        x_init_sss;
        iterations=25,
        ϵ=0.2,
        show_progress=false,
    )
    ms = Matrix(string_sss)
    @test vec(ms[1, :]) ≈ x_init_m[:, 1]
    @test vec(ms[end, :]) ≈ x_init_m[:, end]
end

@testset "String method alg keyword" begin
    Nt_small = 60
    s_small = range(0; stop=1, length=Nt_small)

    xa2 = [-1.0, 0.0]
    xb2 = [1.0, 0.0]
    x_init_m =
        xa2 .* (1 .- s_small)' .+ xb2 .* s_small' .+ [0.0, 0.25] .* sinpi.(s_small)'

    b_nl(x) = [-x[1] + 0.2 * x[2]^3, -0.5 * x[2] - 0.1 * x[1]^3]

    string_default = string_method(b_nl, x_init_m; iterations=20, ϵ=0.3, show_progress=false)
    string_euler = string_method(
        b_nl,
        x_init_m;
        iterations=20,
        ϵ=0.3,
        alg=CriticalTransitions.Euler(),
        show_progress=false,
    )

    @test vec(Matrix(string_default)) ≈ vec(Matrix(string_euler))

    string_tsit5 = string_method(
        b_nl,
        x_init_m;
        iterations=20,
        ϵ=0.3,
        alg=Tsit5(),
        show_progress=false,
    )
    @test norm(vec(Matrix(string_default)) - vec(Matrix(string_tsit5))) > 1e-10
end

@testset "ExtendedPhaseSpace supports alg" begin
    D = 2
    Nt_small = 40
    s_small = range(0; stop=1, length=Nt_small)
    xa2 = [-1.0, 0.0]
    xb2 = [1.0, 0.0]
    x_init_m = xa2 .* (1 .- s_small)' .+ xb2 .* s_small' .+ [0.0, 0.2] .* sinpi.(s_small)'

    H_x_m(x, p) = zeros(size(x))
    H_p_m(x, p) = -x
    sys_m = ExtendedPhaseSpace{false,2}(H_x_m, H_p_m)

    H_x_sss(x, p) = StateSpaceSet(zeros(length(x), D))
    function H_p_sss(x, p)
        m = Matrix(x)
        # Ensure we build a StateSpaceSet from a Nt×D matrix (points in rows),
        # regardless of whether Matrix(StateSpaceSet) returns D×Nt or Nt×D.
        return size(m, 1) == length(x) ? StateSpaceSet(-m) : StateSpaceSet(-permutedims(m))
    end
    sys_sss = ExtendedPhaseSpace{false,2}(H_x_sss, H_p_sss)

    string_euler_m = string_method(
        sys_m,
        x_init_m;
        iterations=15,
        ϵ=0.25,
        alg=CriticalTransitions.Euler(),
        show_progress=false,
    )
    string_tsit5_m = string_method(
        sys_m,
        x_init_m;
        iterations=15,
        ϵ=0.25,
        alg=Tsit5(),
        show_progress=false,
    )

    me = Matrix(string_euler_m)
    mt = Matrix(string_tsit5_m)
    @test vec(me[1, :]) ≈ x_init_m[:, 1]
    @test vec(me[end, :]) ≈ x_init_m[:, end]
    @test vec(mt[1, :]) ≈ x_init_m[:, 1]
    @test vec(mt[end, :]) ≈ x_init_m[:, end]
    @test norm(vec(me) - vec(mt)) > 1e-10

    x_init_sss = StateSpaceSet(x_init_m')
    string_tsit5_sss = string_method(
        sys_sss,
        x_init_sss;
        iterations=15,
        ϵ=0.25,
        alg=Tsit5(),
        show_progress=false,
    )
    ms = Matrix(string_tsit5_sss)
    @test vec(ms[1, :]) ≈ x_init_m[:, 1]
    @test vec(ms[end, :]) ≈ x_init_m[:, end]
end
