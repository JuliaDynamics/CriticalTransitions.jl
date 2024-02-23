using CriticalTransitions, Test

function meier_stein(u, p, t) # out-of-place
    x, y = u
    dx = x-x^3 -10*x*y^2
    dy = -(1+x^2)*y
    [dx, dy]
end
σ = 0.25
sys = StochSystem(meier_stein, [], zeros(2), σ, idfunc, nothing, I(2), "WhiteGauss")

# initial path: parabola
xx = range(-1.0,1.0, length=30)
yy = 0.3 .* (- xx .^ 2 .+ 1)
init = Matrix([xx yy]')

x_i = init[:,1]
x_f = init[:,end]

gm = geometric_min_action_method(sys, init, maxiter=100)
gm = geometric_min_action_method(sys, x_i, x_f, maxiter=100)
path = gm[1][end]
action_val = gm[2][end]
@test all(isapprox.(path[2,:][end-5:end], 0, atol=1e-3))
