using Dierckx, LinearAlgebra, DataInterpolations, Interpolations
using BenchmarkTools
# using CriticalTransitions: anorm, fix_ends, geometric_action

# https://github.com/SciML/DataInterpolations.jl/issues/43
begin # Interpolations.jl is faster for equally spaced ranges
    N = 400
    u = rand(N)
    t = range(0.0, 1.0; length=N)

    x = sort(rand(100)) # sorting required for Interpolations
    println("cubic spline with DataInterpolations")
    interp = DataInterpolations.CubicSpline(u, t)
    @btime $interp($x)
    @btime $interp(0.5)

    println("cubic spline with Interpolations")
    interp2 = Interpolations.CubicSplineInterpolation(t, u)
    @btime $interp2.($x)
    @btime $interp2(0.5)

    println("linear with DataInterpolations")
    interplin = DataInterpolations.LinearInterpolation(u, t)
    @btime $interplin($x)
    @btime $interplin(0.5)

    println("linear with Interpolations")
    interplin2 = Interpolations.LinearInterpolation(t, u)
    @btime $interplin2.($x)
    @btime $interplin2(0.5)
end

function interpolate_Dierckx!(path, α, s)
    α[2:end] .= vec(sqrt.(sum(diff(path; dims=2) .^ 2; dims=1)))
    α .= cumsum(α; dims=1)
    α .= α ./ α[end]
    interp = ParametricSpline(α, path)
    path .= Matrix(interp(s))
    return nothing
end

function interpolate_Interpol!(path, α, s)
    α[2:end] .= vec(sqrt.(sum(diff(path; dims=2) .^ 2; dims=1)))
    α .= cumsum(α; dims=1)
    α .= α ./ α[end]
    path[1, :] .= Interpolations.LinearInterpolation(α, @view path[1, :])(s)
    path[2, :] .= Interpolations.LinearInterpolation(α, @view path[2, :])(s)
    return nothing
end

begin # Interpolations is x3 faster than Dierckx
    using BenchmarkTools
    N = 101
    arclength = 1.0
    xx = range(-1.0, 1.0; length=N)
    yy = 0.3 .* (-xx .^ 2 .+ 1)

    initD = Matrix([xx yy]')
    initI = Matrix([xx yy]')

    alphaD = zeros(N)
    alphaI = zeros(N)
    arc = range(0, 1.0; length=N)

    # interpolate_Dierckx!(initD, alphaD, arc)
    # initD
    # interpolate_Interpol!(initI, alphaI, arc)
    # initI
    # initD .- initI

    @benchmark $interpolate_Dierckx!($initD, $alphaD, $arc)
    @benchmark $interpolate_Interpol!($initI, $alphaI, $arc)
    # Interpolations.LinearInterpolation(t, u)
end
