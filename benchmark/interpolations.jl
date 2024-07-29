using Dierckx, LinearAlgebra
using CriticalTransitions: anorm, fix_ends, geometric_action
using Optim
N = 101
arclength = 1.0
xx = range(-1.0, 1.0; length=N)
yy = 0.3 .* (-xx .^ 2 .+ 1)
init = Matrix([xx yy]')

function interpolate_path(path, N)
    s = zeros(N)
    @inbounds for j in 2:N
        @fastmath s[j] = s[j - 1] + norm(path[:, j] - path[:, j - 1])
    end
    s_length = s / s[end]
    interp = ParametricSpline(s_length, path; k=3)
    return reduce(hcat, [interp(x) for x in range(0, 1.0; length=N)])
end


function interpolate_path!(path, s)
    α[2:end] .= vec(sqrt.(sum(diff(init; dims=2) .^ 2; dims=1)))
    α .= cumsum(α; dims=1)
    α .= α ./ α[end]
    interp = ParametricSpline(α, path)
    path .= Matrix(interp(s))
    return nothing
end

using BenchmarkTools
alpha = zeros(N)
arc = range(0, 1.0; length=N)
@benchmark interpolate_path($init, $N, $arclength)
@benchmark interpolate_path!($alpha, $init, $arc)

path = init
x_i = init[:, 1]
x_f = init[:, end]
N = length(init[1, :])

S(x) = geometric_action(sys, fix_ends(x, x_i, x_f), arclength)
paths = [path]
action = [S(path)]

alpha = zeros(N)
arc = range(0, 1.0; length=N)

update = Optim.optimize(S, path, LBFGS(), Optim.Options(; iterations=1))
path .= Optim.minimizer(update)
Optim.converged(update)
