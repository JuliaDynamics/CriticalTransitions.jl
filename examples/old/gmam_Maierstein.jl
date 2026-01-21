using CriticalTransitions, StaticArrays, LinearAlgebra

using CairoMakie
using MathTeXEngine: get_font
font = (;
    regular=get_font(:regular),
    bold=get_font(:bold),
    italic=get_font(:italic),
    bold_italic=get_font(:bolditalic),
)
using Colors, ColorSchemes

function meier_stein!(du, u, p, t) # in-place
    x, y = u
    du[1] = x - x^3 - 10 * x * y^2
    return du[2] = -(1 + x^2) * y
end
function meier_stein(u, p, t) # out-of-place
    x, y = u
    dx = x - x^3 - 10 * x * y^2
    dy = -(1 + x^2) * y
    return SA[dx, dy]
end
σ = 0.25
sys = CoupledSDEs(meier_stein, idfunc, zeros(2), (), σ)

xx = range(-1.0, 1.0; length=100)
yy = 0.3 .* (-xx .^ 2 .+ 1)
init = Matrix([xx yy]')

gm = geometric_min_action_method(sys, init; maxiters=1000, Stol=1e-6, iter_per_batch=1)

begin
    u_min = -1.1
    u_max = 1.1
    v_min = -0.4
    v_max = 0.4
    res = 100
    u_range = range(u_min, u_max; length=res)
    v_range = range(v_min, v_max; length=res)

    du(u, v) = u - u^3 - 10 * u * v^2
    dv(u, v) = -(1 + u^2) * v
    odeSol(u, v) = Point2f(du(u, v), dv(u, v))

    z = [norm([du(x, y), dv(x, y)]) for x in u_range, y in v_range]
    zmin, zmax = minimum(z), maximum(z)

    stab = [false, true, true]
    fp = [SA[0.0, 0.0], SA[1.0, 0.0], SA[-1.0, 0.0]]
end

begin
    fig = Figure(; size=(600, 400), fontsize=13)
    ax = Axis(
        fig[1, 1];
        xlabel="u",
        ylabel="v",
        aspect=1.4,
        xgridcolor=:transparent,
        ygridcolor=:transparent,
        ylabelrotation=0,
    )

    hm = heatmap!(
        ax, u_range, v_range, z; colormap=:Blues, colorrange=(zmin, zmax), colorscale=sqrt
    )
    Colorbar(fig[1, 2], hm; label="", width=15, ticksize=15, tickalign=1)
    streamplot!(
        ax,
        odeSol,
        (u_min, u_max),
        (v_min, v_max);
        gridsize=(20, 20),
        arrow_size=10,
        stepsize=0.01,
        colormap=[:black, :black],
    )
    colgap!(fig.layout, 7)
    limits!(u_min, u_max, v_min, v_max)
    fig

    [
        scatter!(ax, Point(fp[i]); color=stab[i] > 0 ? :red : :dodgerblue, markersize=10)
        for i in eachindex(fp)
    ]

    lines!(ax, gm[1][1]; linewidth=3, color=:black, linestyle=:dash)
    lines!(ax, gm[1][end]; linewidth=3, color=:orange)
    fig
end
