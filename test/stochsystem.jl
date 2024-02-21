using CriticalTransitions, StaticArrays

function meier_stein(u, p, t) # out-of-place
    x, y = u
    dx = x-x^3 -10*x*y^2
    dy = -(1+x^2)*y
    SA[dx, dy]
end
σ = 0.25

sys = StochSystem(meier_stein, [], zeros(2))

tspan = (0.0, 100.0)
prob = ODEProblem{false}(meier_stein, SA[0.0, 0.0], tspan, ())
kpo = CoupledODEs(prob) # DynamicalSystems method
