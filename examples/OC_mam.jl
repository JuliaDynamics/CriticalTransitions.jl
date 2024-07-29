# using InfiniteOpt, Ipopt, Plots

# T = 50
# opt = Ipopt.Optimizer   # desired solver
# ns = 501;             # number of points in the time grid
# model = InfiniteModel(optimizer_with_attributes(opt));

# @infinite_parameter(model, t in [0, T], num_supports = ns)

# xx(t) = (-1*(1-t/T) + 1*t/T)
# xxx = xx.(supports(t))
# yy(t) = 0.3 .* (- xx(t) .^ 2 .+ 1)
# yyy = yy.(supports(t))
# scatter(xxx, yyy)

# @variable(model, u, Infinite(t), start = xx)
# @variable(model, v, Infinite(t), start = yy)
# du = u - u^3 - 10*u*v^2
# dv = -(1 - u^2)*v

# @objective(model, Min, ∫((∂(u, t) - du)^2 + (∂(v, t) - dv)^2, t))

# @constraint(model, u(0) == -1)
# @constraint(model, v(0) == 0)
# @constraint(model, u(T) == 1)
# @constraint(model, v(T) == 0)

# # @constraint(model, ∂(u, t) == u - u^3 - 10*u*v^2)
# # @constraint(model, ∂(v, t) == -(1 - u^2)*v )

# optimize!(model)

# u_opt = value.(u)
# v_opt = value.(v)
# plot(u_opt, v_opt)
# xlims!(-1.1, 1.1)
# ylims!(-0.1,0.5)
T = 1
using OptimalControl
using NLPModelsIpopt

@def ocp begin
    t ∈ [0, T], time
    u ∈ R, control
    v ∈ R, control
    u(0) == -1
    v(0) == 0
    u(T) == 1
    v(T) == 0
    du = u - u^3 - 10 * u * v^2
    dv = -(1 - u^2) * v
    ∫((∂ₜ(u) - du)^2 + (∂ₜ(v) - dv)^2) → min
end

ocp = Model()                                   # empty optimal control problem
time!(ocp; t0=0, tf=1)                          # initial and final times
state!(ocp, 0)                                  # dimension of the state
control!(ocp, 2)                                # dimension of the control

constraint!(ocp, :initial; lb=[-1, 0], ub=[-1, 0]) # initial condition
constraint!(ocp, :final; lb=[1, 0], ub=[1, 0]) # final condition

dynamics!(ocp, (x, u) -> [x[2], u])           # dynamics of the double integrator

objective!(ocp, :lagrange, (x, u) -> 0.5u^2)    # cost in Lagrange form
control
