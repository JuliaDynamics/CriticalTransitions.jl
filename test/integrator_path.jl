using OrdinaryDiffEq, DynamicalSystemsBase

# Linear ODE which starts at 0.5 and solves from t=0.0 to t=1.0
function b(u::Matrix, p, t)
    1.01*u
end
initu = ones(2,100)
prob = ODEProblem{false}(b, initu, (0.0, Inf))
ds = CoupledODEs(prob)
ds.integ
step!(ds.integ, 0.1)
ds.integ
