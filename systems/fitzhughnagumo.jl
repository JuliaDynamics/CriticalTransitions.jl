"""
Dynamical systems specification file
"""

# fitzhugh_nagumo model
"""
    fitzhugh_nagumo!(du, u, p, t)
In-place definition of the FitzHugh-Nagumo system.

See also [`fitzhugh_nagumo`](@ref).
"""
function fitzhugh_nagumo!(du, u, p, t)
    x, y = u
    ϵ, β, α, γ, κ, I = p

    du[1] = (-α * x^3 + γ * x - κ * y + I) / ϵ
    du[2] = -β * y + x
    return nothing
end

"""
    fitzhugh_nagumo(u, p, t)
Out-of-place definition of the FitzHugh-Nagumo system.

See also [`fitzhugh_nagumo!`](@ref).
"""
function fitzhugh_nagumo(u, p, t)
    x, y = u
    ϵ, β, α, γ, κ, I = p

    dx = (-α * x^3 + γ * x - κ * y + I) / ϵ
    dy = -β * y + x

    return SA[dx, dy]
end

# For backwards compatibility
FitzHughNagumo(u, p, t) = fitzhugh_nagumo(u, p, t)
FitzHughNagumo!(u, p, t) = fitzhugh_nagumo!(u, p, t)

# """
#     fhn_ϵσ(ϵ,σ)
# A shortcut command for returning a CoupledSDEs of the FitzHugh Nagumo system in a default setup with additive isotropic noise.

# This setup fixes the parameters β = 3, α =  γ = κ = 1, I = 0 and leaves the value of the time-scale parameter ϵ as a function argument. The prescribed noise process is additive and isotropic: the variables are peturbed by independently drawn identical Gaussian white noise realisations, with standard deviation σ (the other function argument).
# """
# function fhn_ϵσ(ϵ, σ) # a convenient two-parameter version of the FitzHugh Nagumo system
#     # defining the CoupledSDEs
#     f(u,p,t) = fitzhugh_nagumo(u,p,t);
#     β = 3; α = γ = κ = 1; I = 0; # standard parameters without ϵ (time-scale separation parameter)
#     pf_wo_ϵ = [β, α, γ, κ, I]; # parameter vector without ϵ
#     dim = 2;
#     g = idfunc;
#     pg = nothing;
#     Σ = [1 0; 0 1];
#     process = "WhiteGauss";
#     CoupledSDEs(f, vcat([ϵ], pf_wo_ϵ), dim, σ, g, pg, Σ, process)
# end;

# function fhn_ϵσ_backward(ϵ, σ) # a convenient two-parameter version of the FitzHugh Nagumo system
#     # defining the CoupledSDEs
#     f(u,p,t) = -fitzhugh_nagumo(u,p,t);
#     β = 3; α = γ = κ = 1; I = 0; # standard parameters without ϵ (time-scale separation parameter)
#     pf_wo_ϵ = [β, α, γ, κ, I]; # parameter vector without ϵ
#     dim = 2;
#     g = idfunc;
#     pg = nothing;
#     Σ = [1 0; 0 1];
#     process = "WhiteGauss";
#     CoupledSDEs(f, vcat([ϵ], pf_wo_ϵ), dim, σ, g, pg, Σ, process)
# end;
