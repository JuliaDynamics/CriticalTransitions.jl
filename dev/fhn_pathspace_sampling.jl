"""
Pathspace Langevin MCMC sampling for the FitzHugh-Nagumo model
Author: Reyk
"""

#include("../src/StochSystem.jl")
#include("../src/trajectories/simulation.jl")

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

"""
    fhn_pathspace_sampling(ϵ, σ, T, dz; kwargs...)
Integrates the `FitzHughNagumoSPDE` in virtual time to generate pathspace Langevin MCMC samples.

Unless specified otherwise, the initial condition is a straight line between the fixed points
R and L (for the default FHN parameters).

## Arguments
* `ϵ` (float): time scale parameter
* `σ` (float): noise strength
* `T` (float): total physical time
* `dz` (float): physical time step
* `init=nothing`: initial path. If nothing, `init` is the straight line between the attractors.

## Keyword arguments
* `dt` (float): virtual time step
* `tmax` (float): total virtual time
* `a` (float): first diagonal element of covariance matrix
* `start` (string): Starting point of the path. Either "R" or "L" (default "R")
* `save_every` (int): save the path at every `save_every`-th instance in virtual time
* `β, α, γ, κ, Ι`: FHN parameters

## Returns
A list of two matrices, each with columns = coordinates and rows = virtual time instances.
The first matrix contains the u coordinates, the second matrix the v coordinates.
"""
function fhn_pathspace_sampling(ϵ, σ, T, dz, init=nothing;
    dt=0.001,
    tmax=10.0,
    a=1.0,
    β=3.0,
    α=1.0,
    γ=1.0,
    κ=1.0,
    Ι=0.0,
    save_every=100,
    divterm=true,
    burnin=10)

    # number of points to discretize path
    N = Int(T/dz) + 1
    
    # no noise on boundary conditions
    cov = Matrix(1I(2N))
    cov[1,1] = 0
    cov[N,N] = 0
    cov[N+1,N+1] = 0
    cov[2N, 2N] = 0

    # initial path
    if init == nothing
        # fixed points
        R = [sqrt(2/3), sqrt(2/27)]
        L = - [sqrt(2/3), sqrt(2/27)]

        rr, ll = -range(L[1], R[1], length=N), -range(L[2], R[2], length=N)
        init = vcat(rr,ll)
    end

    # Setup SPDE problem
    spde(eps, sigma) = StochSystem(FitzHughNagumoSPDE, [eps, β, α, γ, κ, Ι, a, Int(divterm)*sigma, dz],
    2N, sqrt(2/dz)*sigma, idfunc, nothing, cov, "WhiteGauss")

    sim = simulate(spde(ϵ,σ), init, dt=dt, tmax=tmax)
    println("$(N) path points, $(round(Int, tmax/dt)) virtual time points, simulation status: $(sim.retcode)")

    sim[1:N,1:save_every:end], sim[N+1:2N,1:save_every:end]
end