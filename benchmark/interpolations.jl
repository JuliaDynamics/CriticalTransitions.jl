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


function interpolate_path!(path, α, s)
    α[2:end] .= vec(sqrt.(sum(diff(path; dims=2) .^ 2; dims=1)))
    α .= cumsum(α; dims=1)
    α .= α ./ α[end]
    interp = ParametricSpline(α, path)
    path .= Matrix(interp(s))
    return nothing
end

using BenchmarkTools
alpha = zeros(N)
arc = range(0, 1.0; length=N)
@benchmark interpolate_path($init, $N)
@benchmark interpolate_path!($init, $alpha, $arc)
