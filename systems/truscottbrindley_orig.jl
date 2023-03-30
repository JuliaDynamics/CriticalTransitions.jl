"""
Dynamical systems specification file
"""

# modified Truscott-Brindley system

"""
    originaltruscottbrindley!(du, u, p, t)
In-place definition of the original Truscott-Brindley system. 

See also [`originaltruscottbrindley`](@ref).
"""
function originaltruscottbrindley!(du, u, p, t)
    P, Z = u
    r, K, Rₘ, α, γ, μ = p[1]

    du[1] = r*P*(1-P/K)-Rₘ*Z*P^2/(α^2+P^2);
    du[2] = γ*Rₘ*Z*P^2/(α^2+P^2)-μ*Z;
end

"""
    originaltruscottbrindley(u, p, t)
Out-of-place definition of the original Truscott-Brindley system. 

See also [`originaltruscottbrindley!`](@ref).
"""
function originaltruscottbrindley(u,p,t)
    P, Z = u
    r, K, Rₘ, α, γ, μ = p[1]

    dP = r*P*(1-P/K)-Rₘ*Z*P^2/(α^2+P^2);
    dZ = γ*Rₘ*Z*P^2/(α^2+P^2)-μ*Z;

    SA[dP, dZ]
end

"""
    rampedoriginaltruscottbrindley!(du, u, p, t)
In-place definition of the ramped original Truscott-Brindley system. 

See also [`rampedoriginaltruscottbrindley`](@ref).
"""
function rampedoriginaltruscottbrindley!(du, u, p, t)

    P, Z, r = u
    K, Rₘ, α, γ, μ, v, Ttrans, Tramp = p[1]

    du[1] = r*P*(1-P/K)-Rₘ*Z*P^2/(α^2+P^2);
    du[2] = γ*Rₘ*Z*P^2/(α^2+P^2)-μ*Z;
    du[3] = t ∈ Ttrans..(Ttrans+Tramp) ? v : 0;

end

"""
    rampedoriginaltruscottbrindley(u, p, t)
Out-of-place definition of the ramped original Truscott-Brindley system. 

See also [`rampedoriginaltruscottbrindley!`](@ref).
"""
function rampedoriginaltruscottbrindley(u, p, t)

    P, Z, r = u
    K, Rₘ, α, γ, μ, v, Ttrans, Tramp = p[1]

    dP = r*P*(1-P/K)-Rₘ*Z*P^2/(α^2+P^2);
    dZ = γ*Rₘ*Z*P^2/(α^2+P^2)-μ*Z;
    dr = t ∈ Ttrans..(Ttrans+Tramp) ? v : 0;

    SA[dP,dZ,dr]
end

"""
    origtb_rσ(r, σ)
A shortcut command for returning a StochSystem of the original Truscott-Brindley system in a default setup with multiplicative anisotropic noise. 
    
This setup fixes the parameters K = 108, Rₘ = 0.7, α = 5.7, γ = 0.05, μ = 0.012 and leaves the value of the growth rate parameter r as a function argument. The prescribed noise process is multiplicative and anisotropic: the first variable is peturbed by Gaussian white noise realisations that are multiplied by the variable's current value; the second variable has no stochastic component. The noise strength σ is left as the remaining function argument.
"""
function origtb_rσ(r, σ) # a convenient three-parameter version of the modifiedtruscottbrindley system 
    f(u,p,t) = originaltruscottbrindley(u,p,t);
    K = 108; Rₘ = 0.7; α = 5.7; γ = 0.05; μ = 0.012; # standard parameters without α (growth rate) and ξ (time-scale separation)
    pf_wo_r = [K, Rₘ, α, γ, μ]; # parameters vector without α or ξ
    dim = 2;
    g(u,p,t) = multiplicative_idx(u,p,t,[true,false]);
    pg = nothing; 
    Σ = [1 0; 0 0];
    process = "WhiteGauss";
    StochSystem(f, vcat([r], pf_wo_r), dim, σ, g, pg, Σ, process)
end;

"""
    rorigtb_vTtrTraσ(v, Ttrans, Tramp, σ)
A shortcut command for returning a StochSystem of the ramped original Truscott-Brindley system in a default setup with multiplicative anisotropic noise. 
    
This setup fixes the parameters K = 108, Rₘ = 0.7, α = 5.7, γ = 0.05, μ = 0.012 and leaves the evolution of the growth rate parameter r as a function argument. The prescribed noise process is multiplicative and anisotropic: the variables are peturbed by Gaussian white noise realisations that are multiplied by their current values - the diffusion matrix is [1 0 0; 0 √γ 0; 0 0 0]. The noise strength σ is left as the remaining function argument.
"""
function rorigtb_vTtrTraσ(v, Ttrans, Tramp, σ) # a convenient three-parameter version of the modifiedtruscottbrindley system 
    f(u,p,t) = rampedoriginaltruscottbrindley(u,p,t);
    K = 108; Rₘ = 0.7; α = 5.7; γ = 0.05; μ = 0.012; # standard parameters without α (growth rate) and ξ (time-scale separation)
    pf_wo_r = [K, Rₘ, α, γ, μ]; # parameters vector without α or ξ
    dim = 3;
    g(u,p,t) = multiplicative_idx(u,p,t,[true,true, false]);
    pg = nothing; 
    Σ = [1 0; 0 √γ 0; 0 0 0];
    process = "WhiteGauss";
    StochSystem(f, vcat(pf_wo_r, [v, Ttrans, Tramp]), dim, σ, g, pg, Σ, process)
end;

