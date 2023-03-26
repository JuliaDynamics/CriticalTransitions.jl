include("planeofbox.jl")

function toattractors(V::Dataset)
    Dict(i => Dataset([V[i]]) for i in 1:length(V))
end

function AVPthreads(X, Y, cds, attractors, p)
    #AVP =  AVPtype(method, p)(cds, attractors);
    ϵ, Ttr, Δt, horizon_limit, mx_chk_lost, diffeq = p; 
    AVP = AttractorsViaProximity(cds,attractors,ϵ;Ttr,Δt,horizon_limit,mx_chk_lost,diffeq)
    A = zeros(length(Y), length(X));
    @Threads.threads for ii ∈ 1:length(X)
        Z = [[X[ii], y] for y ∈ Y];
        A[:,ii] = AVP.(Z)
    end
    A
end

function AVPtype(method::String, p::Vector)
    if method == "defaultCT"
        # (using Ryan's default values for AttractorsViaProximity)
        opt1(x, y) = AttractorsViaProximity(x, y, 0.00005; Ttr = 1000, Δt = 0.001, mx_chk_lost = 100000, diffeq = (alg = Vern9(), abstol = 1e-16, reltol = 1e-16))
    elseif method == "custom"
        # # custom callAVP (using custom values for AttractorsViaProximity)
        # the following seven lines of code would show you how to only specify some arguments, and leave others as preprescribed "standard values", but for now i don't this it's neccesary to implement
        # q = []; # default values of ϵ, trans, dt, horlim, maxitstop, solver
        # for ii ∈ 1:length(p)
        #     if p[ii] != ~
        #         q[ii] = p[ii]
        #     end
        # end
        # ϵ, trans, dt, horlim, maxitstop, solver = q
        ϵ, trans, dt, horlim, maxitstop, solver = p
        opt2(x, y) = AttractorsViaProximity(x, y, ϵ; Ttr = trans, Δt = dt, horizon_limit = horlim, mx_chk_lost = maxitstop, diffeq = solver)
    end
end

"""
    basins(sys::StochSystem, A, B, C, H; kwargs...)
Computes the basins of attraction of StochSystem `sys` on a plane spanned by the points `A`, `B`, `C` and limited by the box `H`. Uses the AttractorsViaProximity function from DynamicalSystems.jl to compute the basins of attraction.

`A`, `B`, `C` are elements of ``\\mathbb{R}^d`` (where ``d`` is the dimension of the  `sys`) and `H` is a hyperrectangle in ``\\mathbb{R}^d``.

The plane is given by ``P_{U,V}\\coloneqq\\{A+u(B-A)+v(C-A): u \\in U,\\, v\\in V\\}`` for some closed and bounded real intervals ``U`` and ``V`` which are selected such that both i) ``P_{U,\\,V} \\subseteq H`` and ii) ``U\\times V\\subseteq\\mathbb{R}^2`` has maximal area, i.e. ``P_{U,\\,V}`` is the "largest" possible plane contained within H. This plane is determined behind the scenes.    
""" 
function basins(sys::StochSystem, A, B, C, H; bstep::Vector = [0.01, 0.01], pstep::Vector = [0.1, 0.1], AVPparams = [0.00005, 1000, 0.001, 1e3, 100000, (alg = Vern9(), abstol = 1e-16, reltol = 1e-16)])

    # here H is a hyperrectangle contained in R^d
    # bstep is the incremental steps you will take across your basin of attraction defined on some plane
    # pstep is a vector of steps (increments) for mechanism behind finding a plane from a box

    projs = plane(A, B, C, H; step = pstep);

    u = projs[1]; v = projs[2];

    # note that the plane returns tuples u, v that tell us the min/max amount to project along the B-A and C-A direction

    # the two dimensional range we loop over 
    X = range(u[1], u[2]; step = bstep[1]); 
    Y = range(v[1], v[2]; step = bstep[2]);

    # the dictionary of attractors

    fps = CriticalTransitions.fixedpoints(sys, H);
    attractors = toattractors(fps[1][fps[3]]);

    h = AVPthreads(X, Y, tocds(sys), attractors, AVPparams);

    [X, Y, attractors, h]

end