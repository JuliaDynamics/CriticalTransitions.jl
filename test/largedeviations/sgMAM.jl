using CriticalTransitions
using ModelingToolkit
using Test

@testset "SgmamSystem KPO" begin
    λ = 3 / 1.21 * 2 / 295
    ω0 = 1.000
    ω = 1.000
    γ = 1 / 295
    η = 0
    α = -1

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


@independent_variables t
D = Differential(t)
sts = @variables u(t) v(t)

eqs = [
    D(u) ~ fu(u, v),
    D(v) ~ fv(u, v)
]
@named sysMTK = System(eqs, t)
sysMTK = structural_simplify(sysMTK)
prob = ODEProblem(sysMTK, sts .=> zeros(2), (0.0, 100.0), (); jac=true)
ds = CoupledODEs(prob)

sys = SgmamSystem{false, 2}(H_x, H_p)
sys′ = SgmamSystem(ds);

Nt = 500  # number of discrete time steps
p_r = rand(2, Nt)
x_r = rand(2, Nt)

@test sys′.H_x(x_r, p_r) ≈ sys.H_x(x_r, p_r)
@test sys′.H_p(x_r, p_r) ≈ sys.H_p(x_r, p_r)

end
