include("planeofbox.jl")

function toattractors(V::Dataset)
    Dict(i => Dataset([V[i]]) for i in 1:length(V))
end

function AVPthreads(X, Y, cds, attractors, method, p)
    callAVP =  AVPtype(method, p)(cds, attractors);
    A = zeros(length(Y), length(X));
    @showprogress for ii ∈ 1:length(X)
        Z = [[X[ii], y] for y ∈ Y];
        A[:,ii] = callAVP.(Z)
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
Computes the basins of attraction of StochSystem `sys` on a plane spanned by the points `A`, `B`, `C` and limited by the box `H`.
"""
function basins(sys::StochSystem, A, B, C, H; fp = nothing, bstep::Vector = [0.01, 0.01], pstep::Vector = [0.1, 0.1], method::String = "defaultCT", AVPparams = [])

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
    sfps = fps[1][fps[3]];
    attractors = toattractors(sfps);

    h = AVPthreads(X, Y, tocds(sys), attractors, method, AVPparams);

    [X, Y, h]

end