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

sys = SgmamSystem(H_x, H_p)

function KPO(x,p,t)
    u, v = x
    return [fu(u, v), fv(u, v)]
end
ds = CoupledSDEs(KPO, zeros(2), ())

# setup
Nt = 500  # number of discrete time steps
s = collect(range(0; stop=1, length=Nt))

xa = [-0.02086931342925046,
0.09908886921365058]
xb = -xa
xsaddle = [0.0, 0.0]

# Initial trajectory
xx = @. (xb[1] - xa[1]) * s + xa[1] + 4 * s * (1 - s) * xsaddle[1]
yy = @. (xb[2] - xa[2]) * s + xa[2] + 4 * s * (1 - s) * xsaddle[2] + 0.01 * sin(2π * s)
x_initial = Matrix([xx yy]')

string = string_method(
    sys, x_initial; iterations=10_000, ϵ=0.5, show_progress=true
)

plot(x_initial[1, :], x_initial[2, :]; label="init", lw=3, c=:black)
plot!(string[1, :], string[2, :]; label="string", lw=3, c=:blue)

string = string_method(
    ds, x_initial; iterations=10_000, ϵ=0.5, show_progress=true
)

plot(x_initial[1, :], x_initial[2, :]; label="init", lw=3, c=:black)
plot!(string[1, :], string[2, :]; label="string", lw=3, c=:blue)

@btime $string_method($sys, $x_initial, iterations=100, ϵ=0.5, show_progress=false) # 7.397 ms (5350 allocations: 12.71 MiB)
@btime $string_method($ds, $x_initial, iterations=100, ϵ=0.5, show_progress=false) # 16.277 ms (253050 allocations: 27.01 MiB)
# @profview string_method(sys, x_initial, iterations=100, ϵ=0.5, show_progress=false)

# The bottleneck is atm at the LinearSolve call to update the x in the new iteration. So the more improve, one needs to write it own LU factorization.
