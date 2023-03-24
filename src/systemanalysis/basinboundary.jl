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
    basinboundary(boa)
Computes the basin boundary for given output `boa` of the [`basins`](@ref) function.

To be further documented.        
"""
basinboundary(boa) = basinboundary(boa[1], boa[2], boa[3])