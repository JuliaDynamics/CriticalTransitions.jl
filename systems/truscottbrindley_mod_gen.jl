function truscottbrindley_mod_gen_det(z,p,t)
    x,y = z; # the system variables
    α,β,ε = p; # the system parameters
    # the drift field
    dx = α*x*(1-x)-x^2*y/(β^2+x^2);
    dy = ε*(x^2*y/(β^2+x^2)-y^2);
    return SVector{2}([dx,dy])
end;