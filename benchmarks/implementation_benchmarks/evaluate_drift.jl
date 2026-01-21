using CriticalTransitions
using Plots
using BenchmarkTools

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

function KPO_SA(x, p, t)
    u, v = x
    return SA[fu(u, v), fv(u, v)]
end
function KPO(x, p, t)
    u, v = x
    return [fu(u, v), fv(u, v)]
end
ds = CoupledSDEs(KPO, zeros(2), ())
ds_sa = CoupledSDEs(KPO_SA, zeros(2), ())

using ModelingToolkit
@independent_variables t
D = Differential(t)
sts = @variables u(t) v(t)

eqs = [D(u) ~ fu(u, v), D(v) ~ fv(u, v)]
@mtkcompile sys1 = System(eqs, t)
prob = ODEProblem(sys1, sts .=> zeros(2), (0.0, 100.0), (); jac=true)
ds = CoupledODEs(prob)
jac = jacobian(ds)
jac([1, 1], (), 0.0)

sgSys′ = SgmamSystem(ds);

p_r = rand(2, Nt)
sgSys′.H_x(x_init_m, p_r) ≈ sys_m.H_x(x_init_m, p_r)
sgSys′.H_p(x_init_m, p_r) ≈ sys_m.H_p(x_init_m, p_r)

@btime sgSys′.H_x($x_init_m, $p_r) # 118.600 μs (5001 allocations: 281.38 KiB)
@btime sys_m.H_x($x_init_m, $p_r) # 5.333 μs (4 allocations: 24.00 KiB)
@btime sgSys′.H_p($x_init_m, $p_r) # 66.200 μs (2001 allocations: 164.19 KiB)
@btime sys_m.H_p($x_init_m, $p_r) # 5.083 μs (4 allocations: 24.00 KiB)
