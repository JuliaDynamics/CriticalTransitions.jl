using InfiniteOpt, Ipopt, Plots

# we have to apply the constraint that T=||d/||B(p)||

T = 9.33
opt = Ipopt.Optimizer   # desired solver
ns = 301;             # number of points in the time grid
model = InfiniteModel(optimizer_with_attributes(opt));

@infinite_parameter(model, t in [0, 1], num_supports = ns)

xx(t) = (-1 * (1 - t) + 1 * t)
xxx = xx.(supports(t))
yy(t) = 0.3 .* (-xx(t) .^ 2 .+ 1)
yyy = yy.(supports(t))
scatter(xxx, yyy)

@variable(model, u, Infinite(t), start = xx)
@variable(model, v, Infinite(t), start = yy)
du = u - u^3 - 10 * u * v^2
dv = -(1 - u^2) * v

@objective(model, Min, T*∫((∂(u, t)/T - du)^2 + (∂(v, t)/T - dv)^2, t)/2)

@constraint(model, u(0) == -1)
@constraint(model, v(0) == 0)
@constraint(model, u(1) == 1)
@constraint(model, v(1) == 0)

@constraint(model, (∂(u, t)^2 + ∂(v, t)^2)/T^2 == du^2 + dv^2)
optimize!(model)

u_opt = value.(u)
v_opt = value.(v)
plot(u_opt, v_opt)
xlims!(-1.1, 1.1)
ylims!(-0.1, 0.5)
