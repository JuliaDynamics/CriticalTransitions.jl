using CriticalTransitions: StochSystem
using CriticalTransitions: simulate, idfunc
using Plots, StaticArrays, LinearAlgebra

############################################################################################
############################################################################################

eps = 0.1 # time scale parameter
sigma = 0.1 # noise strength
T = 15 # total physical time
dz = T/256 # physical time step
dt = 1e-3 # virtual time step
tmax = 20.0 # total virtual time

save_every=100
burnin=10
plot=true

############################################################################################
############################################################################################

a=1.0
β=3.0
α=1.0
γ=1.0
κ=1.0
Ι=0.0

# Central finite difference, first derivative
function central(f, idx, dz)
    (f[idx+1] - f[idx-1])/(2*dz)
end

# Central finite difference, second derivative
function central2(f, idx, dz)
    (f[idx+1] - 2f[idx] + f[idx-1])/(dz^2)
end

"""
    FitzHughNagumoSPDE(u, p, t)
System definition for the pathspace SPDE corresponding to the `fitzhugh_nagumo` system.

The parameter vector must contain the following (in this order):
* ϵ (float): time scale parameter
* β (float): FHN parameter
* α (float): FHN parameter
* γ (float): FHN parameter
* κ (float): FHN parameter
* Ι (float): FHN parameter
* a (float): First diagonal element of covariance matrix 
* σ (float): noise strength
* dz (float): physical time step
"""
function FitzHughNagumoSPDE(u, p, t)
    ϵ, β, α, γ, κ, Ι, a, σ, dz = p[1]

    # N points from 1 to N
    N = Int(length(u)/2)
    du = zeros(length(u))

    for i in 2:N-1
        
        # variable u
        du[i] = (1/a * central2(u,i,dz) + (1+κ/(ϵ*a)) * central(u,i+N,dz)
                    - 1/(ϵ^2*a)*(-α*u[i]^3+γ*u[i]-κ*u[i+N]+Ι)*(-3*α*u[i]^2+γ)
                    + 1*(β*u[i+N]-u[i]) + σ^2*(3*α*u[i]/ϵ))

        # variable v
        du[i+N] = (central2(u,i+N,dz) - (1+κ/(ϵ*a))*central(u,i,dz)
                    + 1*κ/(ϵ^2*a)*(-α*u[i]^3+γ*u[i]-κ*u[i+N]+Ι)
                    + 1*β*(-β*u[i+N]+u[i]))
    end

    SVector{length(u)}(du)
end

# number of points to discretize path
N = Int(T/dz) + 1

cov = Matrix(1I(2N))
cov[1,1] = 0
cov[N,N] = 0
cov[N+1,N+1] = 0
cov[2N, 2N] = 0

# initial path
R = [sqrt(2/3), sqrt(2/27)]
L = - [sqrt(2/3), sqrt(2/27)]
rr, ll = -range(L[1], R[1], length=N), -range(L[2], R[2], length=N)
initpath = vcat(rr,ll)

# Setup SPDE problem
spde(eps, sigma) = StochSystem(FitzHughNagumoSPDE, [eps, β, α, γ, κ, Ι, a, sigma, dz],
2N, sqrt(2/dz)*sigma, idfunc, nothing, cov, "WhiteGauss")


res_u, res_v = [], []
Nstep = round(Int, tmax/dt/save_every)
scatter([R[1],L[1],0], [R[2],L[2],0], xlim=(-1.2,1.2), ylim=(-0.7,0.7), legend=false)
for i in 1:Nstep
    println(i)
    push!(res_u, initpath[1:N])
    push!(res_v, initpath[N+1:2N])
    sim = simulate(spde(eps,sigma), initpath, dt=dt, tmax=tmax/Nstep)
    initpath = sim[:,end]
    
    
    plot!(initpath[1:N], initpath[N+1:2N], markershape=:xcross, color="red")
end

histogram2d(vcat(res_u), vcat(res_v), bins=50)

res_u
