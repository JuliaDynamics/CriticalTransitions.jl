"""
Pathspace Langevin MCMC sampling for the FitzHugh-Nagumo model
Author: Reyk
"""

include("../src/StochSystem.jl")
include("../src/trajectories/simulation.jl")

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
System definition for the pathspace SPDE corresponding to the `FitzHughNagumo` system.

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
        du[i] = 1/a^2 * central2(u,i,dz) + (1+κ/(ϵ*a^2)) * central(u,i+N,dz) - 2/(ϵ^2*a^2)*(-α*u[i]^3+γ*u[i]-κ*u[i+N]+Ι)*(-3*α*u[i]^2+γ) + 2*(β*u[i+N]-u[i]) + σ^2*(3*α*u[i]/ϵ)

        # variable v
        du[i+N] = central2(u,i+N,dz) - (1+κ/(ϵ*a^2))*central(u,i,dz) + 2*κ/(ϵ^2*a^2)*(-α*u[i]^3+γ*u[i]-κ*u[i+N]+Ι) + 2*β*(-β*u[i+N]+u[i])
    end

    SVector{length(u)}(du)
end

"""
    fhn_pathspace_sampling(ϵ, σ, T, dz; kwargs...)
Integrates the `FitzHughNagumoSPDE` in virtual time to generate pathspace Langevin MCMC samples.

The initial condition is a straight line between the fixed points R and L (for the default FHN parameters).

## Arguments
* ϵ (float): time scale parameter
* σ (float): noise strength
* T (float): total physical time
* dz (float): physical time step

## Keyword arguments
* dt (float): virtual time step
* tmax (float): total virtual time
* a (float): first diagonal element of covariance matrix
* start (string): Starting point of the path. Either "R" or "L" (default "R")
* output_every (int): save the path at every `output_every`-th instance in virtual time
* β, α, γ, κ, Ι: FHN parameters

## Returns
A list of two matrices, each with columns = coordinates and rows = virtual time instances.
The first matrix contains the u coordinates, the second matrix the v coordinates.
"""
function fhn_pathspace_sampling(ϵ, σ, T, dz;
    dt=0.001,
    tmax=10.0,
    a=1.0,
    β=3.0,
    α=1.0,
    γ=1.0,
    κ=1.0,
    Ι=0.0,
    start="R",
    output_every=100)

    # number of points to discretize path
    N = Int(T/dz) + 1
    
    # no noise on boundary conditions
    cov = Matrix(1I(2N))
    cov[1,1] = 0
    cov[N,N] = 0
    cov[N+1,N+1] = 0
    cov[2N, 2N] = 0

    # Setup SPDE problem
    spde(eps, sigma) = StochSystem(FitzHughNagumoSPDE, [eps, β, α, γ, κ, Ι, a, sigma, dz], 2N, sqrt(2)*sigma, idfunc, nothing, cov, "WhiteGauss")

    # fixed points
    R = [sqrt(2/3), sqrt(2/27)]
    L = - [sqrt(2/3), sqrt(2/27)]
    
    if start=="R"
        x_i, x_f = R, L
    elseif start=="L"
        x_i, x_f = L, R
    end

    # initial path
    init_path = zeros(2, N)
    for i in 1:N
        init_path[:,i] = x_i .+ (i-1)/(N-1).*(x_f .- x_i)
    end
    init = Vector(vec(init_path'))

    # solve SPDE with Euler-Mayurama
    sim = simulate(spde(ϵ,σ), init, dt=dt, tmax=tmax)
    println("$(N) path points, $(round(Int, tmax/dt)) virtual time points, simulation status: $(sim.retcode)")

    sim[1:N,1:output_every:end], sim[N+1:2N,1:output_every:end]
end