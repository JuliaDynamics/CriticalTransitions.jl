using Optimization,
    OptimizationOptimJL, OptimizationOptimisers, OptimizationOptimJL, OptimizationNLopt  #, Zygote

rosenbrock(x, p) = (p[1] - x[1])^2 + p[2] * (x[2] - x[1]^2)^2
x0 = zeros(2)
p = [1.0, 100.0]

optf = OptimizationFunction(rosenbrock, AutoFiniteDiff())
prob = Optimization.OptimizationProblem(optf, x0, p)
function callback(state, loss)
    state.u .= 0
    return false
end
sol = solve(prob, Optimization.LBFGS()) # changes state
sol = solve(prob, Optimization.LBFGS(); callback=callback) # changes state
sol = solve(prob, Optimisers.Adam(); callback=callback, maxiters=100_000) # changes state
sol = solve(prob, Optimisers.AdamW(); callback=callback, maxiters=100_000) # changes state
sol = solve(prob, NLopt.LD_LBFGS(); callback=callback, maxiters=100_000) # changes state
sol = solve(prob, Optim.LBFGS(); callback=callback, maxiters=100_000) # does not change state
sol = solve(prob, Optim.GradientDescent(); callback=callback, maxiters=100_000) # does not change state
sol = solve(prob, OptimizationOptimJL.LBFGS(); callback=callback) # does not change state
