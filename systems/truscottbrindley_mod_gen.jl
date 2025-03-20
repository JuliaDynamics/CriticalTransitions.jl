function modTBgen_det(z,p,t)
    x,y = z; # the system variables
    α,β,ε = p; # the system parameters
    # the drift field
    dx = α*x*(1-x)-x^2*y/(β^2+x^2);
    dy = ε*(x^2*y/(β^2+x^2)-y^2);
    return SVector{2}([dx,dy])
end;

function truscottbrindley_mod_gen_det(p::Vector{Float64},diffeq;
    t0::Float64 = 0.,
    kwargs...)
    # finding the real equilibrium with greatest x-value
    α,β = p[1:2];
    xeq = roots([α*β^4,-α*β^4,2*α*β^2,-1-2*α*β^2,α,-α]);
    xeqreal = real.(xeq[abs.(imag.(xeq).<1.e-12)]);
    xᵣ = findmax(xeqreal)[1];
    zᵣ = [xᵣ,α*(1-xᵣ)*(β^2+xᵣ^2)/xᵣ];
    # returning the CoupledODEs
    return CoupledODEs(modTBgen_det,zᵣ,p;diffeq,t0,kwargs...)
end

