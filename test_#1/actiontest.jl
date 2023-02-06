include("../src/CriticalTransitions.jl")
using .CriticalTransitions

sys = StochSystem(FitzHughNagumo, [0.1,3.,1.,1.,1.,0.], 2, 0.08)

sys.f

path = zeros(2,10)
for i in 1:10
    path[1,i] = sqrt(2/3) + (i-1)/9*(-2*sqrt(2/3))
    path[2,i] = sqrt(2/27) + (i-1)/9*(-2*sqrt(2/27))
end
path

sys.f([0.1,0.2], [sys.pf, sys.pg], 0)

size(path, 2)

using LinearAlgebra

function anorm(vec, A; square=false)
    normsquared = dot(vec, A * vec)
    if square
        return normsquared
    else
        return sqrt(normsquared)
    end
end

anorm(path[:,1], [1 0; 0 1])

dot(path[:,1], path[:,1])

A = inv([1 0; 0 1])

sqrt(dot(path[:,1], A * path[:,1]))

function p(sys::StochSystem)
    [sys.pf, sys.pg]
end

function fw_integrand(sys::StochSystem, path, time)
    v = zeros(size(path))
    v[:,1] .= (path[:,2] .- path[:,1]) ./ (time[2] .- time[1])
    v[:,end] .= (path[:,end] .- path[:,end-1]) ./ (time[end] .- time[end-1])
    for i in 2:(size(path, 2) - 1)
        v[:, i] .= (path[:,i+1] .- path[:,i-1]) ./ (time[i+1] .- time[i-1])
    end
    A = inv(sys.Σ)
    sqnorm = zeros(size(path, 2))
    for i in 1:size(path, 2)
        drift = sys.f(path[:,i], p(sys), time[i])
        sqnorm[i] = anorm(v[:,i] - drift, A, square=true)
    end
    sqnorm
end;


function fw_action(sys::StochSystem, path, time)
    integrand = fw_integrand(sys, path, time)
    S = 0
    for i in 1:(size(path, 2) - 1)
        S += (integrand[i+1] + integrand[i])/2 * (time[i+1]-time[i])
    end
    S/2
end;


time = range(0,15,266)
fw_action(sys, path, time)

using HDF5

file = h5open("research/software/CriticalTransitions.jl/test/maps_langevinmcmc_eps1.0e-01.h5", "r")
path1 = file["path_saddle"][:,:]'
path2 = file["path_detached"][:,:]'
time1 = range(0,15,266)

fw_action(sys, path2, time1)

using ForwardDiff

function om_action(sys::StochSystem, path, time)
    integrand = fw_integrand(sys, path, time)
    S = 0
    for i in 1:(size(path, 2) - 1)
        S += (integrand[i+1] + integrand[i])/2 * (time[i+1]-time[i])
        S += sys.σ^2 * (div_b(sys, path[:,i+1]) + div_b(sys, path[:,i]))/2 * (time[i+1]-time[i])
    end
    S/2
end

function div_b(sys::StochSystem, x)
    b(x) = sys.f(x, p(sys), 0)
    tr(ForwardDiff.jacobian(b, x))
end

div_b(sys, path[:,2])


b(x) = sys.f(x, p(sys), 0)
ForwardDiff.jacobian(b, path[:,end])

# through saddle
fw_action(sys, path1, time1)
om_action(sys, path1, time1)

# detached
fw_action(sys, path2, time1)
om_action(sys, path2, time1)

# straight line
fw_action(sys, path, time)
om_action(sys, path, time)
