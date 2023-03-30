"""
Dynamical systems specification file
"""

# FitzHughNagumo model

"""
    FitzHughNagumo!(du, u, p, t)
In-place definition of the FitzHugh-Nagumo system.

See also [`FitzHughNagumo`](@ref).
"""
function FitzHughNagumo!(du, u, p, t)
    x, y = u
    ϵ, β, α, γ, κ, I = p[1]

    du[1] = (-α*x^3 + γ*x - κ*y + I)/ϵ
    du[2] = -β*y + x
end

"""
    FitzHughNagumo(u, p, t)
Out-of-place definition of the FitzHugh-Nagumo system.

See also [`FitzHughNagumo!`](@ref).
"""
function FitzHughNagumo(u,p,t)
    x, y = u
    ϵ, β, α, γ, κ, I = p[1]

    dx = (-α*x^3 + γ*x - κ*y + I)/ϵ
    dy = -β*y + x

    SA[dx, dy]
end


"""
    fhn_ϵσ(ϵ,σ)
A shortcut command for returning a StochSystem of the FitzHugh Nagumo system in a default setup with additive isotropic noise. 
    
This setup fixes the parameters β = 3, α =  γ = κ = 1, I = 0 and leaves the value of the time-scale parameter ϵ as a function argument. The prescribed noise process is additive and isotropic: the variables are peturbed by independently drawn identical Gaussian white noise realisations, with standard deviation σ (the other function argument).
"""
function fhn_ϵσ(ϵ, σ) # a convenient two-parameter version of the FitzHugh Nagumo system 
    # defining the StochSystem
    f(u,p,t) = FitzHughNagumo(u,p,t);
    β = 3; α = γ = κ = 1; I = 0; # standard parameters without ϵ (time-scale separation parameter)
    pf_wo_ϵ = [β, α, γ, κ, I]; # parameter vector without ϵ
    dim = 2;
    g = idfunc;
    pg = nothing; 
    Σ = [1 0; 0 1];
    process = "WhiteGauss";
    StochSystem(f, vcat([ϵ], pf_wo_ϵ), dim, σ, g, pg, Σ, process)
end;

function fhn_ϵσ_backward(ϵ, σ) # a convenient two-parameter version of the FitzHugh Nagumo system 
    # defining the StochSystem
    f(u,p,t) = -FitzHughNagumo(u,p,t);
    β = 3; α = γ = κ = 1; I = 0; # standard parameters without ϵ (time-scale separation parameter)
    pf_wo_ϵ = [β, α, γ, κ, I]; # parameter vector without ϵ
    dim = 2;
    g = idfunc;
    pg = nothing; 
    Σ = [1 0; 0 1];
    process = "WhiteGauss";
    StochSystem(f, vcat([ϵ], pf_wo_ϵ), dim, σ, g, pg, Σ, process)
end;