include("planeofbox.jl")
include("../StochSystem.jl")
include("stability.jl")

function toattractors(V::Dataset)
    Dict(i => Dataset([V[i]]) for i in 1:length(V))
end

function AVPthreads(X, Y, cds, attractors, callAVP)
    # the following loop ensures that X is the vector with the greatest length
    if length(Y)>length(X)
        Y, X = [X, Y];
    end 
    A = zeros(length(Y), length(X));
    @Threads.threads for ii ∈ 1:length(X)
        Z = [[X[ii], y] for y ∈ Y];
        A[:,ii] = callAVP(cds, attractors).(Z)
    end
    A
end

function callAVP(method::String, p::Vector)
    if method == "defaultCT"
        # (using Ryan's default values for AttractorsViaProximity)
        foo(cds, attractors) = AttractorsViaProximity(cds, attractors, 0.00005; Ttr = 1000, Δt = 0.001, mx_chk_lost = 100000, diffeq = (alg = Vern9(), abstol = 1e-16, reltol = 1e-16));
    elseif method == "defaultDS"
        # (using George's default values for AttractorsViaProximity)
        ϵ = p[1];
        foo(cds, attractors) = AttractorsViaProximity(cds, attractors, ϵ);
    elseif method == "custom"
        # # custom callAVP (using custom values for AttractorsViaProximity)
        ϵ, trans, dt, horlim, mxitstop, solver = p
        foo(cds, attractors) = AttractorsViaProximity(cds, attractors, ϵ; Ttr = trans, Δt = dt, horizon_limit = horlim, mx_chk_lost = maxitstp, diffeq = solver)
    end
    foo
end

function basins(sys::StochSystem, A, B, C, H; fp = nothing, stp = [0.01, 0.01], method::String = "defaultCT", AVPparams = [])

    # here H is a hyperrectangle contained in R^d

    projs = plane(A, B, C, H);

    u = projs[1]; v = projs[2];

    # note that the plane returns tuples u, v that tell us the min/max amount to project along the B-A and C-A direction

    # the two dimensional range we loop over 
    X = range(u[1], u[2]; step = stp[1]); 
    Y = range(v[1], v[2]; step = stp[2]);

    # the dictionary of attractors

    fps = CriticalTransitions.fixedpoints(sys, H);
    sfps = fps[1][fps[3]];
    attractors = toattractors(sfps);

    h = AVPthreads(X, Y, tocds(sys), attractors, callAVP(method, AVPparams))

end