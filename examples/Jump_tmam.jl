using InfiniteOpt, Ipopt, Plots

# The basic idea of tMAM is yo remove the optimization parameter T by replacing it with an optimal linear time scaling
# Define the submodel
opt = Ipopt.Optimizer   # desired solver
ns = 201;             # number of points in the time grid
begin
    submodel = InfiniteModel(optimizer_with_attributes(opt, "tol"=>1e-9))
    # Setup up the normalized time parameter τ = t / tf
    @infinite_parameter(submodel, τ in [0, 1], num_supports = ns)
    @variable(submodel, 0.1 <= tf <= 100, start = 1)

    xx(t) = (-1 * (1 - t) + 1 * t)
    yy(t) = 0.5 .* (-xx(t) .^ 2 .+ 1)

    u = parameter_function(xx, τ)
    v = parameter_function(yy, τ)
    du = u - u^3 - 10 * u * v^2
    dv = -(1 - u^2) * v

    @objective(submodel, Min, tf * ∫((∂(u, τ) / tf - du)^2 + (∂(v, τ) / tf - dv)^2, τ) / 2)

    optimize!(submodel)
    tf_prime = value(tf)
end

model = InfiniteModel(optimizer_with_attributes(opt));

@infinite_parameter(model, τ in [0, 1], num_supports = ns)

xxx = xx.(supports(τ))
yyy = yy.(supports(τ))
scatter(xxx, yyy)

@variable(model, u, Infinite(τ), start = xx)
@variable(model, v, Infinite(τ), start = yy)
du = u - u^3 - 10 * u * v^2
dv = -(1 - u^2) * v

@objective(
    model,
    Min,
    tf_prime * ∫((∂(u, τ) / tf_prime - du)^2 + (∂(v, τ) / tf_prime - dv)^2, τ) / 2
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
