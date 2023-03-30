include("planeofbox.jl")

toattractors(V::Dataset) = Dict(i => StateSpaceSet([V[i]]) for i in 1:length(V))

"""
    basins(sys::StochSystem, A, B, C, H; kwargs...)
Computes the basins of attraction of StochSystem `sys` on a plane spanned by the distinct points `A`, `B`, `C` and limited by the box `H`. Uses the [`AttractorsViaProximity`](https://juliadynamics.github.io/Attractors.jl/v1.2/attractors/#Attractors.AttractorsViaProximity) function from [`DynamicalSystems.jl`](https://juliadynamics.github.io/DynamicalSystems.jl/stable/) to compute the basins of attraction.

`A`, `B`, `C` are elements of ``\\mathbb{R}^d`` (where ``d`` is the dimension of the  `sys`) and `H` is a hyperrectangle in ``\\mathbb{R}^d``.

The plane is given by ``P_{U,V} := \\{A+u(B-A)+v(C-A)\\in\\mathbb{R}^d: u \\in U,\\, v\\in V\\}`` for some closed and bounded real intervals ``U`` and ``V`` which are selected such that both 
* ``P_{U,\\,V} \\subseteq H``, and 
* ``U\\times V\\subseteq\\mathbb{R}^2`` has maximal area, 
i.e. ``P_{U,\\,V}`` is the "largest" possible plane contained within H. This plane is determined behind the scenes.    

This function returns a four-dimensional vector. The first two entries are discretised versions of the interval ``U`` and ``V`` (as defined above, of lengths ``\\ell_U,\\,\\ell_V`` respectively); the third entry is a dictionary of the attractors (stable equilibria) of the system within `H`, and the final entry is an ``\\ell_V\\times\\ell_U`` matrix of integers that group the initial conditions (written in terms of ``A+u(B-A)+v(C-A)`` where ``u\\in U`` and ``v\\in V``) by which attractor they will in time converge to. 

## Keyword arguments 
* `bstep = [0.01, 0.01]`: a vector of length two whose elements respectively specify the length of the incremental steps taken across each dimension in the discretisation of your plane
* `pstep = [0.1, 0.1]`: a vector of length two whose elements give the increments of the mesh that the maximisation process of finding a plane from a box is taken over (for more information see the source code of the function `plane` in the `src/systemanalysis/planeofbox.jl` file)
* `ϵ_mapper = 0.01`: `ϵ` parameter of [`AttractorsViaProximity`](https://juliadynamics.github.io/Attractors.jl/v1.2/attractors/#Attractors.AttractorsViaProximity)
* `kwargs...`: keyword arguments passed to the [`AttractorsViaProximity`](https://juliadynamics.github.io/Attractors.jl/v1.2/attractors/#Attractors.AttractorsViaProximity) function (namely, `Ttr, Δt, horizon_limit, mx_chk_lost`)
""" 
function basins(sys::StochSystem, A, B, C, H; bstep::Vector = [0.01, 0.01], pstep::Vector = [0.1, 0.1], ϵ_mapper=0.01, kwargs...)

    U, V = plane(A, B, C, H; step = pstep); # the intervals for the maximal plane P_{U,V}

    X = range(U[1], U[2]; step = bstep[1]); # the range of projections along the B-A direction 
    Y = range(V[1], V[2]; step = bstep[2]); # the range of projections along the C-A direction
    h = zeros(length(Y),length(X)); # an array to store the attractors for each initial condition

    # defining the parameters required to run the AttractorsViaProximity function
    fps = CriticalTransitions.fixedpoints(sys, H);
    attractors = toattractors(fps[1][fps[3]]);

    # running the AttractorsViaProximity function using parallel computing
    @Threads.threads for ii ∈ tqdm(1:length(X))
        Z = [A+X[ii]*(B-A)+y*(C-A) for y ∈ Y]; # a row of initial conditions
        h[:,ii] = AttractorsViaProximity(CoupledODEs(sys), attractors, ϵ_mapper; kwargs...).(Z)
    end

    [X, Y, attractors, h]

end