using InfiniteOpt, Ipopt, Plots

T = 50
opt = Ipopt.Optimizer   # desired solver
ns = 501;             # number of points in the time grid
model = InfiniteModel(optimizer_with_attributes(opt));

@infinite_parameter(model, t in [0, T], num_supports = ns)

xx(t) = (-1*(1-t/T) + 1*t/T)
xxx = xx.(supports(t))
yy(t) = 0.3 .* (- xx(t) .^ 2 .+ 1)
yyy = yy.(supports(t))
scatter(xxx, yyy)

@variable(model, u, Infinite(t), start = xx)
@variable(model, v, Infinite(t), start = yy)
du = u - u^3 - 10*u*v^2
dv = -(1 - u^2)*v

@objective(model, Min, ∫((∂(u, t) - du)^2 + (∂(v, t) - dv)^2, t))

@constraint(model, u(0) == -1)
@constraint(model, v(0) == 0)
@constraint(model, u(T) == 1)
@constraint(model, v(T) == 0)

# @constraint(model, ∂(u, t) == u - u^3 - 10*u*v^2)
# @constraint(model, ∂(v, t) == -(1 - u^2)*v )

optimize!(model)

u_opt = value.(u)
v_opt = value.(v)
plot(u_opt, v_opt)
xlims!(-1.1, 1.1)
ylims!(-0.1,0.5)
