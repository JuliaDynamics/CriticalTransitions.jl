"""
Dynamical systems specification file
"""

# modified Truscott-Brindley model for Ryan's mad tipping explorations!!!

function modifiedtruscottbrindley!(du, u, p, t)
    P, Z = u
    α, β, γ, P₁, Z₁, ξ = p[1]

    du[1] = P₁*(α*(P/P₁)*(1-β*(P/P₁))-γ*(Z/Z₁)*(P/P₁)^2/(1+(P/P₁)^2));
    du[2] = ξ*Z₁* ((Z/Z₁)*(P/P₁)^2/(1+(P/P₁)^2)-(Z/Z₁)^2);
end

function modifiedtruscottbrindley(u,p,t)
    P, Z = u
    α, β, γ, P₁, Z₁, ξ = p[1]

    dP = P₁*(α*(P/P₁)*(1-β*(P/P₁))-γ*(Z/Z₁)*(P/P₁)^2/(1+(P/P₁)^2));
    dZ = ξ*Z₁* ((Z/Z₁)*(P/P₁)^2/(1+(P/P₁)^2)-(Z/Z₁)^2);

    SVector{2}([dP, dZ])
end

function modtb_αξσ(α, ξ, σ) # a convenient three-parameter version of the modifiedtruscottbrindley system 
    # defining the StochSystem
    f(u,p,t) = modifiedtruscottbrindley(u,p,t);
    β = 5/112; γ = 112/(45*0.0525); P₁ = β; Z₁ = 5/6; # standard parameters without α (growth rate) and ξ (time-scale separation)
    pf_wo_αξ = [β, γ, P₁, Z₁];
    dim = 2;
    #σ = 0.1; 
    g(u,p,t) = multiplicative_idx(u,p,t,[true,false]);
    pg = nothing; 
    Σ = [1 0; 0 0];
    process = "WhiteGauss";
    StochSystem(f, vcat([α], pf_wo_αξ, [ξ]), dim, σ, g, pg, Σ, process)
end;