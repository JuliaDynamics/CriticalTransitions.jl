"""
Dynamical systems specification file
"""

# modified Truscott-Brindley system

"""
    modifiedtruscottbrindley!(du, u, p, t)
In-place definition of the modified Truscott-Brindley system. 

See also [`modifiedtruscottbrindley`](@ref).
"""
function modifiedtruscottbrindley!(du, u, p, t)
    P, Z = u
    α, β, γ, P₁, Z₁, ξ = p[1]

    du[1] = P₁*(α*(P/P₁)*(1-β*(P/P₁))-γ*(Z/Z₁)*(P/P₁)^2/(1+(P/P₁)^2));
    du[2] = ξ*Z₁* ((Z/Z₁)*(P/P₁)^2/(1+(P/P₁)^2)-(Z/Z₁)^2);
end

"""
    modifiedtruscottbrindley(u, p, t)
Out-of-place definition of the modified Truscott-Brindley system. 

See also [`modifiedtruscottbrindley!`](@ref).
"""
function modifiedtruscottbrindley(u,p,t)
    P, Z = u
    α, β, γ, P₁, Z₁, ξ = p[1]

    dP = P₁*(α*(P/P₁)*(1-β*(P/P₁))-γ*(Z/Z₁)*(P/P₁)^2/(1+(P/P₁)^2));
    dZ = ξ*Z₁* ((Z/Z₁)*(P/P₁)^2/(1+(P/P₁)^2)-(Z/Z₁)^2);

    SVector{2}([dP, dZ])
end

"""
    rampedmodifiedtruscottbrindley!(du, u, p, t)
In-place definition of the ramped modified Truscott-Brindley system. 

See also [`rampedmodifiedtruscottbrindley`](@ref).
"""
function rampedmodifiedtruscottbrindley!(du, u, p, t)

    P, Z, α = u
    β, γ, P₁, Z₁, ξ, v, Ttrans, Tramp = p[1]

    du[1] = P₁*(α*(P/P₁)*(1-β*(P/P₁))-γ*(Z/Z₁)*(P/P₁)^2/(1+(P/P₁)^2));
    du[2] = ξ*Z₁* ((Z/Z₁)*(P/P₁)^2/(1+(P/P₁)^2)-(Z/Z₁)^2);
    du[3] = t ∈ Ttrans..(Ttrans+Tramp) ? v : 0;

end

"""
    rampedmodifiedtruscottbrindley(u, p, t)
Out-of-place definition of the ramped modified Truscott-Brindley system. 

See also [`rampedmodifiedtruscottbrindley!`](@ref).
"""
function rampedmodifiedtruscottbrindley(u, p, t)

    P, Z, α = u
    β, γ, P₁, Z₁, ξ, v, Ttrans, Tramp = p[1]

    dP = P₁*(α*(P/P₁)*(1-β*(P/P₁))-γ*(Z/Z₁)*(P/P₁)^2/(1+(P/P₁)^2));
    dZ = ξ*Z₁* ((Z/Z₁)*(P/P₁)^2/(1+(P/P₁)^2)-(Z/Z₁)^2);
    dα = t ∈ Ttrans..(Ttrans+Tramp) ? v : 0;

    SVector{3}([dP,dZ,dα])

end

"""
    modtb_αξσ(α, ξ, σ)
A shortcut command for returning a StochSystem of the modified Truscott-Brindley system in a default setup with multiplicative anisotropic noise. 
    
This setup fixes the parameters β = 5/112, γ = 112/2.3625, P₁ = β, Z₁ = 5/6 and leaves the values of the parameters α and ξ as function arguments. The prescribed noise process is multiplicative and anisotropic: the first variable is peturbed by Gaussian white noise realisations that are multiplied by the variable's current value; the second variable has no stochastic component. The noise strength σ is left as the remaining function argument.
"""
function modtb_αξσ(α, ξ, σ) # a convenient three-parameter version of the modifiedtruscottbrindley system 
    f(u,p,t) = modifiedtruscottbrindley(u,p,t);
    β = 5/112; γ = 112/(45*0.0525); P₁ = β; Z₁ = 5/6; # standard parameters without α (growth rate) and ξ (time-scale separation)
    pf_wo_αξ = [β, γ, P₁, Z₁]; # parameters vector without α or ξ
    dim = 2;
    g(u,p,t) = multiplicative_idx(u,p,t,[true,false]);
    pg = nothing; 
    Σ = [1 0; 0 0];
    process = "WhiteGauss";
    StochSystem(f, vcat([α], pf_wo_αξ, [ξ]), dim, σ, g, pg, Σ, process)
end;

function modtb_αξσ_backward(α, ξ, σ) # a convenient three-parameter version of the modifiedtruscottbrindley system 
    f(u,p,t) = -modifiedtruscottbrindley(u,p,t);
    β = 5/112; γ = 112/(45*0.0525); P₁ = β; Z₁ = 5/6; # standard parameters without α (growth rate) and ξ (time-scale separation)
    pf_wo_αξ = [β, γ, P₁, Z₁]; # parameters vector without α or ξ
    dim = 2;
    g(u,p,t) = multiplicative_idx(u,p,t,[true,false]);
    pg = nothing; 
    Σ = [1 0; 0 0];
    process = "WhiteGauss";
    StochSystem(f, vcat([α], pf_wo_αξ, [ξ]), dim, σ, g, pg, Σ, process)
end;

"""
    rmodtb_ξvTtrTraσ(ξ, v, Ttrans, Tramp, σ)
A shortcut command for returning a StochSystem of the ramped modified Truscott-Brindley system in a default setup with multiplicative anisotropic noise. 
    
This setup fixes the parameters β = 5/112, γ = 112/2.3625, P₁ = β, Z₁ = 5/6 and leaves the time-scale separation parameter and the evolution of the growth rate parameter r as a function argument. The prescribed noise process is multiplicative and anisotropic: the variables are peturbed by Gaussian white noise realisations that are multiplied by their current values - the diffusion matrix is [1 0 0; 0 0 0; 0 0 0]. The noise strength σ is left as the remaining function argument.
"""
function rmodtb_ξvTtrTraσ(ξ, v, Ttrans, Tramp, σ) # a convenient three-parameter version of the modifiedtruscottbrindley system 
    f(u,p,t) = rampedmodifiedtruscottbrindley(u,p,t);
    β = 5/112; γ = 112/(45*0.0525); P₁ = β; Z₁ = 5/6; # standard parameters without α (growth rate) and ξ (time-scale separation)
    pf_wo = [β, γ, P₁, Z₁]; # parameters vector without α or ξ
    dim = 3;
    g(u,p,t) = multiplicative_idx(u,p,t,[true,false,false]);
    pg = nothing; 
    Σ = [1 0 0; 0 0 0; 0 0 0];
    process = "WhiteGauss";
    StochSystem(f, vcat(pf_wo, [ξ, v, Ttrans, Tramp]), dim, σ, g, pg, Σ, process)
end;
