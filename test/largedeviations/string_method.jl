using CriticalTransitions

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

function H_x_sss(x, p) # ℜ² → ℜ²
    u, v = eachcol(x)
    pu, pv = eachcol(p)

    H_u = @. pu * dfudu(u, v) + pv * dfvdu(u, v)
    H_v = @. pu * dfudv(u, v) + pv * dfvdv(u, v)
    return StateSpaceSet([H_u H_v]) |> vec
end
function H_p_sss(x, p) # ℜ² → ℜ²
    u, v = eachcol(x)
    pu, pv = eachcol(p)

    H_pu = @. pu + fu(u, v)
    H_pv = @. pv + fv(u, v)
    return StateSpaceSet([H_pu H_pv]) |> vec
end

sys_sss = SgmamSystem{false,2}(H_x_sss, H_p_sss)

x_init_sss = StateSpaceSet([xx yy])
sys_sss.H_p(x_init_sss, 0 * Matrix(x_init_sss))

function H_x(x, p) # ℜ² → ℜ²
    u, v = eachrow(x)
    pu, pv = eachrow(p)

    H_u = @. pu * dfudu(u, v) + pv * dfvdu(u, v)
    H_v = @. pu * dfudv(u, v) + pv * dfvdv(u, v)
    return Matrix([H_u H_v]')
end
function H_p(x, p) # ℜ² → ℜ²
    u, v = eachrow(x)
    pu, pv = eachrow(p)

    H_pu = @. pu + fu(u, v)
    H_pv = @. pv + fv(u, v)
    return Matrix([H_pu H_pv]')
end

sys_m = SgmamSystem{false,2}(H_x, H_p)

x_init_m = Matrix([xx yy]')

@testset "Integrators" begin
    using CriticalTransitions: updateEuler!, updateRK4!

    @testset "Euler" begin
        ref1 = deepcopy(x_init_m)
        ref2 = deepcopy(x_init_sss)
        updateEuler!(ref1, sys_m, 0.1)
        updateEuler!(ref2, sys_sss, 0.1)
        @test ref1 ≈ Matrix(Matrix(ref2)')
    end

    @testset "RK4" begin
        ref1 = deepcopy(x_init_m)
        ref2 = deepcopy(x_init_sss)
        updateRK4!(ref1, sys_m, 0.1)
        # updateRK4!(ref2, sys_sss, 0.1) #broken
        @test ref1 ≈ Matrix(Matrix(ref2)')

    end
end

@testset "StateSpaceSet vs Matrix" begin
    string_sss = string_method(
        sys_sss, x_init_sss; iterations=10_000, ϵ=0.5, show_progress=false
    )

    string_m = string_method(sys_m, x_init_m; iterations=10_000, ϵ=0.5, show_progress=false)

    @test vec(string_m) ≈ vec(string_sss)
end
