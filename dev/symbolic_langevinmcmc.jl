#include("../src/StochSystem.jl")

function langevinmcmc_spde(u, p, t)
    
    # here u should be a vector of length N × dim where N is the number of (physical) time points in the discretisation of [0,Tphys]

    ## in the parameters vector we need:
    # the dimension of the system, dim
    # the noise strength of the system, σ
    # the (invertible) covariance matrix, Σ
    # the physical time step, Δz
    # the state variables in symbolic form, x_{i=1}^dim, state_vars
    # the jacobian of the drift, jacobian
    # the penultimate term as a callable function, grad_dot_term
    # the final term as a function, grad_div_term

    # t represents the virtual time

    dim, σ, Σ, Δz, state_vars, jacobian, grad_dot_term, grad_div_term = p[1]
    
    N = Int(length(u)/dim); # the number of components across the path
    du = zeros(length(u)); # a dummy storage for the derivatives

    for ii ∈ 2:N-1
        
        dzxii = [central(u[(jj-1)*N+1:jj*N], ii, Δz) for jj ∈ 1:dim]; # the derivative ∂x/∂z at physical time z = ii*Δz 
        dzzxii = [central2(u[(jj-1)*N+1:jj*N], ii, Δz) for jj ∈ 1:dim]; # the derivative ∂^2x/∂z^2 at physical time z = ii*Δz
        xii = [u[(jj-1)*N+ii] for jj ∈ 1:dim]; # the state variable state_vars at z = ii*Δz

        J = substitute.(jacobian, (Dict(zip(state_vars,xii)),)); # the jacobian evaluated at xii 
        Jxii = Symbolics.value.(J); # making this usable
        f = substitute.(grad_dot_term, (Dict(zip(state_vars,xii)),)); # the ∇⟨b, inv(Σ)*b⟩ term evaluated at xii
        fxii = Symbolics.value.(f); # making this usable
        g = substitute.(grad_div_term, (Dict(zip(state_vars,xii)),)); # the ∇(∇⋅b) term evaluated at xii
        gxii = Symbolics.value.(g); # making this usable

        du[[(jj-1)*N+ii for jj ∈ 1:dim]] = inv(Σ)*dzzxii - (inv(Σ)*Jxii - Jxii'*inv(Σ))*dzxii - (1/2)*fxii - (σ^2/2)*gxii; # the deterministic portion of the SPDE function evaluated at xii

    end

    SVector{length(u)}(du)

end

function symbolise_spde(sys::StochSystem)

    # creating symbolic state variables
    state_vars = @eval @variables $([Symbol("x$i") for i=1:length(sys.u)]...);

    # next, the deterministic drift of the system in symbolic form
    drift = sys.f(state_vars,[sys.pf],@variables t)

    # next, the jacobian of the drift in symbolic form
    jacobian = Symbolics.jacobian(drift,state_vars);

    # next, the penultimate deterministic term in the LangevinMCMC SPDE
    grad_dot_term = Symbolics.gradient(dot(drift, inv(sys.Σ)*drift), state_vars);
    
    # finally, the last deterministic term in the LangevinMCMC SPDE
    grad_div_term = Symbolics.gradient(tr(jacobian), state_vars); # since the divergence is purely the trace of the jacobian
    
    return state_vars, jacobian, grad_dot_term, grad_div_term

end

function jacobian(sys::StochSystem)
    # creating symbolic state variables
    state_vars = @eval @variables $([Symbol("x$i") for i=1:length(sys.u)]...);

    # next, the deterministic drift of the system in symbolic form        
    drift = sys.f(state_vars,[sys.pf],@variables t);
    
    # next, the jacobian of the drift in symbolic form
    jacobian = Symbolics.jacobian(drift,state_vars)

end

function jacobian(sys::StochSystem, x::Vector)

    # loading the Jacobian in symbolic form
    jacobian1 = jacobian(sys);

    # calling the Jacobian at the state value x
    state_vars = @eval @variables $([Symbol("x$i") for i=1:length(sys.u)]...);
    J = substitute.(jacobian1, (Dict(zip(state_vars,x)),)); # the jacobian evaluated at x
    Jx = Symbolics.value.(J);

    return jacobian1, Jx

end
    

"""
    langevinmcmc(sys::StochSystem, init; kwargs...)
Runs Langevin MCMC bridge sampling for a given system `sys` and initial path `init`, using
Symbolics.jl to compute the functional derivative.

Returns a three-dimensional array of size `(M, length(sys.u), N)`, where `M` is the number of
virtual time steps, `N` the number of physical time steps, and `length(sys.u)` the number of system
dimensions.

To be further documented in following versions.

## Keyword arguments
* `T = 15`: physical path time
* `tmax = 250`: virtual time duration
* `Δt = 1e-3`: virtual time step
* `noise = sys.σ`: noise strength of additive spatio-temporal noise term
* `annealing = false`: whether to enable simulated annealing
* `showprogress = true`: whether to print a progress bar
"""
function langevinmcmc(sys::StochSystem, init;    
    T = 15,         # physical end time
    tmax = 250.0,   # virtual end time
    Δt = 1e-3,      # virtual time step
    noise = sys.σ,  # noise strength of additive term
    annealing = false,  # simulated annealing
    showprogress= true)

    N = size(init, 2); # the number of columns in the initial condition gives the number of components of the path 
    Δz = T/(N-1); # the physical time step
    t = range(0, T; length = N); # the physical time range

    Nstep = round(Int64, tmax/Δt); # the number of mcmc iterations

    paths = zeros(Nstep+1, length(sys.u), N); # three-dimensional array to store the paths
    paths[1,:,:] = init; # storing the initial condition
    t = 0; # the current virtual time value 

    # defining the number of functions used in langevin mcmc 
    println("... Symbolizing SPDE")
    state_vars, jacobian, grad_dot_term, grad_div_term = symbolise_spde(sys);

    p = [length(sys.u), sys.σ, sys.Σ, Δz, state_vars, jacobian, grad_dot_term, grad_div_term]; # langevinmcmcSPDE parameters 

    println("... Initializing sampling")

    iterator = showprogress ? tqdm(1:Nstep) : 1:Nstep
    for it ∈ 1:Nstep

        if annealing
            _noise = noise*(Nstep-it)/Nstep
        else
            _noise = noise
        end

        update = langevinmcmc_spde(vec(paths[it,:,:]'), [p], t); # this gives the virtual-time derivatives for all components on the current path

        # we need to split update (dim*N × 1) into a dim × N matrix
        update_mod = zeros(length(sys.u), N)
        for jj ∈ 1:length(sys.u)
            update_mod[jj,:] = update[(jj-1)*N+1:jj*N];
        end 

        paths[it+1,:,:] = paths[it,:,:] .+ (Δt * update_mod .+ _noise * sqrt(2*Δt/Δz)*randn(length(sys.u),N)); # the Euler-Maruyama step 
        
        # reseting the fixed end points
        for jj ∈ 1:length(sys.u)
            paths[it+1,jj,[1,end]] = paths[it,jj,[1,end]]
        end

        t += Δt;

    end

    paths

end

function stochtolangevinmcmcstoch(sys::StochSystem, Tphys::Float64, Δz::Float64) 

    N = 1+ceil(Int64,Tphys/Δz); # number of path points

    f(u,p,t) = langevinmcmc_spde(u,p,t); # defining the f-field of StochSystem
    p = vcat([length(sys.u), sys.σ, sys.Σ, Δz],[symbolise_spde(sys)[i] for i ∈ 1:length(symbolise_spde(sys))]); # defining the p-field of StochSystem
    dim = length(sys.u)*N; # the dim-field of StochSystem
    σ = √(2/Δz)*sys.σ; # defining the σ-field of StochSystem

    z = ones(dim); 
    z[reduce(vcat,[[(ii-1)*N+1, ii*N] for ii ∈ 1:length(sys.u)])] .= zeros(2*length(sys.u)); # this is a modified version of idfunc with zeros at the beginning and end of each length(sys.u) portion
    g(u,p,t) = SVector{dim}(z); # defining the g-field of StochSystem

    pg = [];
    Σ = I(dim);
    process = "WhiteGauss";

    StochSystem(f, p, dim, σ, g, pg, Σ, process)

end
