using InfiniteOpt, Ipopt, Plots

opt = Ipopt.Optimizer   # desired solver
ns = 101;             # number of points in the time grid
model = InfiniteModel(optimizer_with_attributes(opt));

@infinite_parameter(model, τ in [0, 1], num_supports = ns)

xx(t) = (-1 * (1 - t) + 1 * t)
yy(t) = 0.3 .* (-xx(t) .^ 2 .+ 1)

@variable(model, u, Infinite(τ), start = xx)
@variable(model, v, Infinite(τ), start = yy)
@variable(model, T)
du = u - u^3 - 10 * u * v^2
dv = -(1 - u^2) * v
∂u = @deriv(u, τ)
∂v = @deriv(v, τ)


@objective(
    model,
    Min,
    T * ∫((∂u - du)^2 + (∂v - dv)^2, τ) / 2
)

@constraint(model, u(0) == -1)
@constraint(model, v(0) == 0)
@constraint(model, u(1) == 1)
@constraint(model, v(1) == 0)

optimize!(model)

u_opt = value.(u)
v_opt = value.(v)
plot(u_opt, v_opt)
xlims!(-1.1, 1.1)
ylims!(-0.1, 0.5)
