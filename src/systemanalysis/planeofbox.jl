function intmultints(v)
    # here v is a vector of intervals
    w = v[1];
    for ii ∈ 2:length(v)
        w = w ∩ v[ii];
    end
    w
end

function vthreshold(u, A, B, C, box)
    U = B - A;
    V = C - A;
    ϕ₁ = [];
    for ii ∈ 1:length(A)
        if V[ii] ≠ 0
            vals = ((box[ii].lo - A[ii] - u[1]*U[ii])/V[ii], (box[ii].hi - A[ii] - u[1]*U[ii])/V[ii])
            low = minimum(vals);
            high = maximum(vals);
            ϕ₁ = push!(ϕ₁, low..high)
        end
    end
    φ₁ = intmultints(ϕ₁);
    ϕ₂ = [];
    for ii ∈ 1:length(A)
        if V[ii] ≠ 0
            vals = ((box[ii].lo - A[ii] - u[2]*U[ii])/V[ii], (box[ii].hi - A[ii] - u[2]*U[ii])/V[ii])
            low = minimum(vals);
            high = maximum(vals);
            ϕ₂ = push!(ϕ₂, low..high)
        end
    end
    φ₂ = intmultints(ϕ₂);
    φ = intersect(φ₁, φ₂)
    return [φ.lo, φ.hi] 
end

function normalise(A::Vector)
    return A/norm(A)
end

function plane(A, B, C, box; step::Vector = [0.1, 0.1])

    # box is a hyperrectangle in R^d         
    # here A is the starting vector, and B and C are vectors that determine the direction of the plane

    if A ∉ box 
        @warn("Please enter a starting vector that is contained within the prescribed box")
    elseif !all([length(B), length(C), length(box)] .== length(A))
        @warn("Please enter vectors A, B, C and a set 'box' that are of the same dimension")
    elseif A == zeros(length(A)) == B || A == zeros(length(A)) == C
        @warn("Please enter vectors such that neither A=B=0 or A=C=0 hold") 
    elseif normalise(B-A) == normalise(C-A) # i.e. in particular B =/= C
        @warn("Please enter vectors B and C that are not collinear")
    elseif A == B || A == C
        @warn("Please enter distinct vectors A, B, and C.")
    else ## main code

        d = length(A); # the (shared) dimension that we are working within 
        ϕₜ = []; # here i should define an empty vector in a better way! i'm later going to fill it with Interval{Float64}s
            
        @Threads.threads for ii ∈ 1:d
            if A[ii] ≠ B[ii]
                vals = ((box[ii].lo-A[ii])/(B[ii]-A[ii]), (box[ii].hi-A[ii])/(B[ii]-A[ii]));
                low = minimum(vals);
                high = maximum(vals);
                ϕₜ = push!(ϕₜ, low..high); # store the range of permitted t values for A + t (B-A) in the iith component 
            end
        end
        
        φₜ = intmultints(ϕₜ)

        ## so right now we have φₜ and we wish to subsequently loop over these 
        
        U₁ = range(0, φₜ.lo; step = -step[1]);
        U₂ = range(0, φₜ.hi; step = step[2]);

        V = [[] for i ∈ 1:length(U₁), j ∈ 1:length(U₂)]

        @Threads.threads for ii ∈ 1:length(U₁)
            us = [[U₁[ii], u] for u ∈ U₂]; # the candidate for umin and umax
            vs = vthreshold.(us, Ref(A), Ref(B), Ref(C), Ref(box));
            V[ii,:] = [[[(U₁[ii], U₂[jj]), (vs[jj][1], vs[jj][2])], (vs[jj][2] - vs[jj][1])*(U₂[jj] - U₁[ii])] for jj ∈ 1:length(U₂)]
        end

        W = [V[ii,jj][2] for ii ∈ 1:size(V,1), jj ∈ 1:size(V,2)];

        q = findmax(W)[2];

        V[q][1];

    end
end

# function plane(A, B, C, box; step::Vector = [0.1, 0.1])

#     # box is a hyperrectangle in R^d         
#     # here A is the starting vector, and B and C are vectors that determine the direction of the plane

#     if A ∉ box 
#         @warn("Please enter a starting vector that is contained within the prescribed box")
#     elseif !all([length(B), length(C), length(box)] .== length(A))
#         @warn("Please enter vectors A, B, C and a set 'box' that are of the same dimension")
#     elseif A == zeros(length(A)) == B || A == zeros(length(A)) == C
#         @warn("Please enter vectors such that neither A=B=0 or A=C=0 hold") 
#     elseif normalise(B-A) == normalise(C-A)
#         @warn("Please enter vectors B and C that are not collinear")
#     elseif A == B || A == C || B = C
#         @warn("Please enter distinct vectors A, B, and C.")
#     else ## main code
#         d = length(A); # the (shared) dimension that we are working within 
#         ## firstly we deal in the (B-A) direction 
#         ϕₜ = []; # here i should define an empty vector in a better way! i'm later going to fill it with Interval{Float64}s
#         if A ≠ B
#             @Threads.threads for ii ∈ 1:d
#                 if A[ii] ≠ B[ii]
#                     vals = ((box[ii].lo-A[ii])/(B[ii]-A[ii]), (box[ii].hi-A[ii])/(B[ii]-A[ii]));
#                     low = minimum(vals);
#                     high = maximum(vals);
#                     ϕₜ = push!(ϕₜ, low..high); # store the range of permitted t values for A + t (B-A) in the iith component 
#                 end
#             end
#         else # i.e. A = B, this is a separate interpretation!!! for when you are considering A+t*B but B=A
#             @Threads.threads for ii ∈ 1:d
#                 if A[ii] ≠ 0
#                     vals = (box[ii].lo/A[ii], box[ii].hi/A[ii]) 
#                     low = minimum(vals)
#                     high = maximum(vals)
#                     ϕₜ = push!(ϕₜ, low..high);
#                 end
#             end
#         end
#         φₜ = intmultints(ϕₜ)
#         ## next we deal with the (C-A) direction (actually irrelevant for the current way the code is written)
#         ϕₛ = [] # here i should define an empty vector in a better way! i'm later going to fill it with Interval{Float64}s
#         if A ≠ C
#             @Threads.threads for ii ∈ 1:d
#                 if A[ii] ≠ C[ii]
#                     vals = ((box[ii].lo-A[ii])/(C[ii]-A[ii]), (box[ii].hi-A[ii])/(C[ii]-A[ii]));
#                     low = minimum(vals);
#                     high = maximum(vals);
#                     ϕₛ = push!(ϕₛ, low..high); # store the range of permitted t values for A + s (C-A) in the ith component 
#                 end
#             end
#         else # i.e. A = C 
#             @Threads.threads for ii ∈ 1:d
#                 if A[ii] ≠ 0
#                     vals = (box[ii].lo/A[ii], box[ii].hi/A[ii]) 
#                     low = minimum(vals)
#                     high = maximum(vals)
#                     ϕₛ = push!(ϕₛ, low..high);
#                 end
#             end
#         end
#         φₛ = intmultints(ϕₛ)
#         ## so right now we have φₜ, φₛ and we wish to subsequently loop over these 
#         U₁ = range(0, φₜ.lo; step = -step[1]);
#         U₂ = range(0, φₜ.hi; step = step[2]);

#         #println([length(U₁), length(U₂)])

#         ## method one

#         # V = [[] for i ∈ 1:length(U₁), j ∈ 1:length(U₂)] 

#         # @showprogress for ii ∈ 1:length(U₁), jj ∈ 1:length(U₂)
#         #     v = vthreshold([U₁[ii], U₂[jj]], A, B, C, box)
#         #     V[ii,jj] = [[(U₁[ii], U₂[jj]), (v[1], v[2])], (v[2]-v[1])*(U₂[jj]-U₁[ii])]
#         # end

#         ## method two

#         V = [[] for i ∈ 1:length(U₁), j ∈ 1:length(U₂)]

#         @Threads.threads for ii ∈ 1:length(U₁)
#             us = [[U₁[ii], u] for u ∈ U₂]; # the candidate for umin and umax
#             vs = vthreshold.(us, Ref(A), Ref(B), Ref(C), Ref(box));
#             V[ii,:] = [[[(U₁[ii], U₂[jj]), (vs[jj][1], vs[jj][2])], (vs[jj][2] - vs[jj][1])*(U₂[jj] - U₁[ii])] for jj ∈ 1:length(U₂)]
#         end

#         W = [V[ii,jj][2] for ii ∈ 1:size(V,1), jj ∈ 1:size(V,2)];

#         q = findmax(W)[2];

#         V[q][1];

#         ## method three

#         # V = [[] for i ∈ 1:length(U₁)*length(U₂)] 

#         # Us = [[U₁[1], u] for u ∈ U₂];
#         # for jj ∈ 2:length(U₂)
#         #     Us = vcat(Us, [[U₁[jj], u] for u ∈ U₂]);
#         # end

#         # @showprogress for ii ∈ 1:length(Us) 
#         #     v = vthreshold(Us[ii], A, B, C, box);
#         #     V[ii] = [[(Us[ii][1], Us[ii][2]), (v[1], v[2])], (v[2]-v[1])*(Us[ii][2]-Us[ii][1])]
#         # end
        
#         ## method four

#         # Us = [[u₁, u₂] for u₁ ∈ U₁, u₂ ∈ U₂];

#         # Vs = vthreshold.(Us, Ref(A), Ref(B), Ref(C), Ref(box));

#         # V = [[[(Us[ii,jj][1], Us[ii,jj][2]), (Vs[ii,jj][1], Vs[ii,jj][2])], (Vs[ii,jj][2]-Vs[ii,jj][1])*(Us[ii,jj][2]-Us[ii,jj][1])] for ii ∈ 1:length(U₁), jj ∈ 1:length(U₂)];

#         ## implement parallel computing + broadcasting and check whether this is faster!!
#         ## compare with the matrix broadcasting method

#         # corners  = zeros(1,5); 
#         # @showprogress for u₁ ∈ U₁, u₂ ∈ U₂
#         #     for v₁ ∈ V₁, v₂ ∈ V₂
#         #         if [A .+ u₁*B + v₁.*C, A .+ u₁*B + v₂.*C, A .+ u₂*B + v₁.*C, A .+ u₂*B + v₂.*C] ⊆ box
#         #             # the corners of the plane (and hence the plane itself) are contained in the box
#         #             corners = vcat(corners, [u₁, u₂, v₁, v₂, (u₂-u₁)*(v₂-v₁)]')
#         #         end
#         #     end
#         # end
#         # corners = corners[2:end,:]
#     #return V
#     end
# end