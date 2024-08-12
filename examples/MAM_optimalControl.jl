using OptimalControl
using NLPModelsIpopt
using Plots

T = 50

@def action begin
    t ∈ [0, T], time
    x ∈ R², state
    u ∈ R², control
    x(0) == [-1, 0]
    x(T) == [1, 0]
    ẋ(t) == u(t)
    ∫( sum((u(t) - f(x(t))).^2) ) → min
end

f(u, v) = [u - u^3 - 10*u*v^2,  -(1 - u^2)*v]
f(x) = f(x...)

x1(t) = -(1 - t / T) + t / T
x2(t) = 0.3(-x1(t)^2 + 1)
x(t) = [x1(t), x2(t)]
u(t) = f(x(t))
init = (state=x, control=u)

sol = solve(action; init=init, grid_size=50)
sol = solve(action; init=sol, grid_size=1000) # grid refinement

plot(sol)
sol |> propertynames
 @which plot!(sol)
state =  sol.state.(sol.times)
plot(first.(state), last.(state))
