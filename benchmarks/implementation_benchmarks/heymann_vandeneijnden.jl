using CriticalTransitions
using BenchmarkTools

function meier_stein(u, p, t) # out-of-place
    x, y = u
    dx = x - x^3 - 10 * x * y^2
    dy = -(1 + x^2) * y
    return SA[dx, dy]
end

σ = 0.25
sys = CoupledSDEs(meier_stein, zeros(2); noise_strength=σ)

N = 200
xx = range(-1.0, 1.0; length=N)
yy = 0.3 .* (-xx .^ 2 .+ 1)
init = Matrix([xx yy]')

@btime CriticalTransitions.geometric_gradient_step($sys, $init, $N; tau=0.1)
@btime geometric_min_action_method(
    $sys,
    $init;
    maxiters=10,
    show_progress=false,
    optimizer=GeometricGradient(),
    stepsize=0.1,
)
