include("planeofbox.jl")

toattractors(V::Dataset) = Dict(i => Dataset([V[i]]) for i in 1:length(V))

# function AVPthreads(A, B, C, X, Y, cds, attractors, p)
#     # AVP =  AVPtype("defaultCT", p)(cds, attractors);
#     ϵ, Ttr, Δt, horizon_limit, mx_chk_lost, diffeq = p; 
#     #AVP = AttractorsViaProximity(cds,attractors,ϵ;Ttr,Δt,horizon_limit,mx_chk_lost,diffeq)
#     A = zeros(length(Y), length(X));
#     @Threads.threads for ii ∈ tqdm(1:length(X))
#         #Z = [[X[ii], y] for y ∈ Y];
#         Z = [A+X[ii]*(B-A)+y*(C-A) for y ∈ Y]; # a row of initial conditions
#         A[:,ii] = AttractorsViaProximity(cds,attractors,ϵ;Ttr,Δt,horizon_limit,mx_chk_lost,diffeq).(Z)
#     end
#     A
# end

# function AVPtype(method::String, p::Vector)
#     if method == "defaultCT"
#         # (using Ryan's default values for AttractorsViaProximity)
#         opt1(x, y) = AttractorsViaProximity(x, y, 0.00005; Ttr = 1000, Δt = 0.001, mx_chk_lost = 100000, diffeq = (alg = Vern9(), abstol = 1e-16, reltol = 1e-16))
#     elseif method == "custom"
#         # # custom callAVP (using custom values for AttractorsViaProximity)
#         # the following seven lines of code would show you how to only specify some arguments, and leave others as preprescribed "standard values", but for now i don't this it's neccesary to implement
#         # q = []; # default values of ϵ, trans, dt, horlim, maxitstop, solver
#         # for ii ∈ 1:length(p)
#         #     if p[ii] != ~
#         #         q[ii] = p[ii]
#         #     end
#         # end
#         # ϵ, trans, dt, horlim, maxitstop, solver = q
#         ϵ, trans, dt, horlim, maxitstop, solver = p
#         opt2(x, y) = AttractorsViaProximity(x, y, ϵ; Ttr = trans, Δt = dt, horizon_limit = horlim, mx_chk_lost = maxitstop, diffeq = solver)
#     end
# end

"""
    basins(sys::StochSystem, A, B, C, H; kwargs...)
Computes the basins of attraction of StochSystem `sys` on a plane spanned by the distinct points `A`, `B`, `C` and limited by the box `H`. Uses the AttractorsViaProximity function from DynamicalSystems.jl to compute the basins of attraction.

`A`, `B`, `C` are elements of ``\\mathbb{R}^d`` (where ``d`` is the dimension of the  `sys`) and `H` is a hyperrectangle in ``\\mathbb{R}^d``.

The plane is given by ``P_{U,V} := \\{A+u(B-A)+v(C-A)\\in\\mathbb{R}^d: u \\in U,\\, v\\in V\\}`` for some closed and bounded real intervals ``U`` and ``V`` which are selected such that both 
* ``P_{U,\\,V} \\subseteq H``, and 
* ``U\\times V\\subseteq\\mathbb{R}^2`` has maximal area, 
i.e. ``P_{U,\\,V}`` is the "largest" possible plane contained within H. This plane is determined behind the scenes.    

This function returns a four-dimensional vector. The first two entries are discretised versions of the interval ``U`` and ``V`` (as defined above, of lengths ``\\ell_U,\\,\\ell_V`` respectively); the third entry is a dictionary of the attractors (stable equilibria) of the system within `H`, and the final entry is an ``\\ell_V\\times\\ell_U`` matrix of integers that group the initial conditions (written in terms of ``A+u(B-A)+v(C-A)`` where ``u\\in U`` and ``v\\in V``) by which attractor they will in time converge to. 

## Keyword arguments 
* `bstep = [0.01, 0.01]`: a vector of length two whose elements respectively specify the length of the incremental steps taken across each dimension in the discretisation of your plane
* `pstep = [0.1, 0.1]`: a vector of length two whose elements give the increments of the mesh that the maximisation process of finding a plane from a box is taken over (for more information see the source code of the function `plane` in the src/systemanalysis/planeofbox.jl file)
* `AVPparams = [0.00005, 1000, 0.001, 1e3, 100000, (alg = Vern9(), abstol = 1e-16, reltol = 1e-16)]`: a vector of length six specifying parameters for the DynamicalSystems.AttractorsViaProximity function (namely, `ϵ, Ttr, Δt, horizon_limit, mx_chk_lost, maxitstop, diffeq`)
""" 
function basins(sys::StochSystem, A, B, C, H; bstep::Vector = [0.01, 0.01], pstep::Vector = [0.1, 0.1], AVPparams = [0.00005, 1000, 0.001, 1e3, 100000, (alg = Vern9(), abstol = 1e-16, reltol = 1e-16)])

    U, V = plane(A, B, C, H; step = pstep); # the intervals for the maximal plane P_{U,V}

    X = range(U[1], U[2]; step = bstep[1]); # the range of projections along the B-A direction 
    Y = range(V[1], V[2]; step = bstep[2]); # the range of projections along the C-A direction
    h = zeros(length(Y),length(X)); # an array to store the attractors for each initial condition

    # defining the parameters required to run the AttractorsViaProximity function
    fps = CriticalTransitions.fixedpoints(sys, H);
    attractors = toattractors(fps[1][fps[3]]);
    ϵ, Ttr, Δt, horizon_limit, mx_chk_lost, diffeq = AVPparams; 

    # running the AttractorsViaProximity function using parallel computing
    @Threads.threads for ii ∈ tqdm(1:length(X))
        Z = [A+X[ii]*(B-A)+y*(C-A) for y ∈ Y]; # a row of initial conditions
        h[:,ii] = AttractorsViaProximity(tocds(sys),attractors,ϵ;Ttr,Δt,horizon_limit,mx_chk_lost,diffeq).(Z)
    end

    [X, Y, attractors, h]

end