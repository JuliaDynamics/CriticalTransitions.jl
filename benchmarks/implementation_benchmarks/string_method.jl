using CriticalTransitions
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

function sss()
    function H_x(x, p) # ℜ² → ℜ²
        u, v = eachcol(x)
        pu, pv = eachcol(p)

        H_u = @. pu * dfudu(u, v) + pv * dfvdu(u, v)
        H_v = @. pu * dfudv(u, v) + pv * dfvdv(u, v)
        return StateSpaceSet([H_u H_v])
    end
    function H_p(x, p) # ℜ² → ℜ²
        u, v = eachcol(x)
        pu, pv = eachcol(p)

        H_pu = @. pu + fu(u, v)
        H_pv = @. pv + fv(u, v)
        return StateSpaceSet([H_pu H_pv])
    end

    sys_sss = ExtendedHamiltonianSystem(H_x, H_p)

    x_init_sss = StateSpaceSet([xx yy])
    return sys_sss, x_init_sss
end

function m()
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

    sys_m = ExtendedHamiltonianSystem(H_x, H_p)

    x_init_m = Matrix([xx yy]')
    return sys_m, x_init_m
end

sys_sss, x_init_sss = sss()
sys_m, x_init_m = m()

@benchmark string_method($sys_sss, $x_init_sss)
@benchmark string_method($sys_m, $x_init_m)

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

@btime string_method($ds_sa, $x_init_sss) # 77.119 ms (48693 allocations: 86.62 MiB)
@btime string_method($ds, $x_init_sss) # 151.156 ms (1044693 allocations: 162.61 MiB)
@btime string_method($sys_m, $x_init_m) # 150.201 ms (51540 allocations: 127.86 MiB)
@btime string_method($ds_sa, $x_init_m) # 178.669 ms (1039716 allocations: 158.70 MiB)
@btime string_method($sys_sss, $x_init_sss) #  206.504 ms (1087517 allocations: 186.73 MiB)
@btime string_method($ds, $x_init_m) # 244.689 ms (2533716 allocations: 272.68 MiB)
