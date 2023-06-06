function planetophase(A::Vector, B::Vector, C::Vector, proj::Vector)
    A + proj[1]*(B-A) + proj[2]*(C-A)
end

"""
    basinboundary(X, Y, h; kwargs...)
Computes the basin boundary for given output `X, Y, h` of the [`basins`](@ref) function.

To be further documented.
"""
function basinboundary(X, Y, h; coords::String = "plane", A::Vector = [], B::Vector = [], C::Vector = [])

    bb = [];

    for ii ∈ 1:size(h,2)     
        for jj ∈ 1:size(h,1)-1 # traversing through each column of the basins of attraction
           if h[jj,ii] != h[jj+1,ii] # if two points on top of each other are not in the same boa
                if !any(v->v==-1, [h[jj,ii], h[jj+1, ii]]) # each of these points go to detected attractor
                    np = [X[ii], (Y[jj]+Y[jj+1])/2]; # the value
                    bb = push!(bb,np);
                else # i.e. one of the points has value negative 1
                    if h[jj,ii] == -1 # i.e. the saddle point is here
                        np = [X[ii], Y[jj]];
                    else # the saddle point is at the point (X[ii], Y[jj+1])
                        np = [X[ii], Y[jj+1]];
                    end
                    bb = push!(bb, np);
                end
            end
        end
    end

    bb = unique(bb); # to account for the saddles counted twice

    if coords == "plane"
        planecoords = [bb[jj][ii] for ii ∈ 1:2, jj ∈ 1:length(bb)] # convert into nice matrix form
    elseif coords == "phase"
        bb_plane = planetophase.(Ref(A), Ref(B), Ref(C), bb); 
        phasecoords = [bb_plane[jj][ii] for ii ∈ 1:length(A), jj ∈ 1:length(bb)]    
    end
    
end

"""
    basboundary(sys::StochSystem, xrange::Vector, yrange::Vector, xspacing::Float64, attractors::Vector; kwargs...)

This function computes the basin boundary.
"""

function basboundary(sys::StochSystem, xrange::Vector, yrange::Vector, xspacing::Float64, attractors::Vector;
    eps1=1e-4,
    ϵ_mapper=0.1,
    dt_mapper=1.0e-3,
    solver=Vern9(),
    maxit=1e+5)

    #N = convert(Int64,(xrange[2]-xrange[1])/xspacing)
    N = round(Int64,(xrange[2]-xrange[1])/xspacing)
    bby = zeros(N+1)
    xx = range(xrange[1], xrange[2], length=N+1)

    Threads.@threads for i ∈ 1:N+1   
        #println(i)
        u1,u2 = bisect_to_edge2(sys, [xx[i],yrange[1]], [xx[i],yrange[2]], attractors,
        eps1=eps1,
        ϵ_mapper=ϵ_mapper,
        dt_mapper=dt_mapper,
        solver=solver,
        maxit=maxit)

        bby[i] = ((u1 + u2)/2)[2]
    end

    xx,bby

end

"""
    basinboundary(boa)
Computes the basin boundary for given output `boa` of the [`basins`](@ref) function.

To be further documented.        
"""
basinboundary(boa) = basinboundary(boa[1], boa[2], boa[4])