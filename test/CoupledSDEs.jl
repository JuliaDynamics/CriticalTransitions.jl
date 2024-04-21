using CriticalTransitions

function meier_stein(u, p, t) # out-of-place
    x, y = u
    dx = x - x^3 - 10 * x * y^2
    dy = -(1 + x^2) * y
    SA[dx, dy]
end
σ = 0.25

CoupledSDEs(meier_stein, (u,p,t) -> σ*idfunc(u,p,t), zeros(2))
