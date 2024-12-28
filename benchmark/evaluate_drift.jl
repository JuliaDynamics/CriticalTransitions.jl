# System setup - FitzHugh-Nagumo model
p = [1.0, 3.0, 1.0, 1.0, 1.0, 0.0] # Parameters (ϵ, β, α, γ, κ, I)
σ = 0.2 # noise strength
sys = CoupledSDEs(fitzhugh_nagumo, zeros(2), p; noise_strength=σ)

A = inv(covariance_matrix(sys))
T, N = 2.0, 100

x_i = SA[sqrt(2 / 3), sqrt(2 / 27)]
x_f = SA[0.0, 0.0]

path = reduce(hcat, range(x_i, x_f; length=N))
time = range(0.0, T; length=N)

b(x) = drift(sys, x)
b(x[:,2])

using DynamicalSystemsBase
b.(StateSpaceSet(path'))
x = deepcopy(path)
function mod!(x)
    for i in 1:size(x)[2]
        x[:,i] .= b(x[:,i])
    end
end
x = deepcopy(path)
function mod1!(sys, x)
    for i in 1:size(x)[2]
        x[:,i] .= sys.integ.f(x[:,i], sys.p0, 0)
    end
end
function mod2!(sys, x)
    for i in 1:size(x)[2]
        x[:,i] = sys.integ.f(x[:,i], sys.p0, 0)
    end
end
x = deepcopy(path)
@benchmark mod!($x)
x = deepcopy(path)
@benchmark mod1!($sys,$x)
x = deepcopy(path)
@benchmark mod2!($sys,$x)
