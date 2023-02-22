#include("../src/StochSystem.jl")

"""
    edgetracking(sys::StochSystem, u1::State, u2::State, attractors::Vector; kwargs...)
Runs the edge tracking algorithm.

## Input 
* `sys`: dynamical system of type [`StochSystem`](@ref)
* `u1`, `u2`: initial states; must belong to different basins of attraction
* `attractors`: vector of state vectors corresponding to the stable fixed points of `sys`

## Keyword arguments
* `eps1 = 1e-9`: tolerance for bisection distance
* `eps2 = 1e-8`: tolerance for divergence of trajectories before re-bisecting
* `converge = 1e-5`: convergence criterion for M state accuracy (Euclidean distance)
* `dt = 0.01`: integration time step 
* `tmax = Inf`: maximum integration time of parallel trajectories until re-bisection 
* `ϵ_mapper = 0.1`: distance threshold for AttractorMapper
* `solver = Vern9()`: ODE solver from `DifferentialEquations.jl`
* `maxit = 100`: maximum number of iterations before algorithm stops 
* `verbose = true`: print status updates during run
* `output_all = false`: if `false`, returns M state, else returns all points of the track
* `kwargs...`: additional keyword arguments of `AttractorsViaProximity` may be passed

## Returns
If `output_all`, a single state vector corresponding to the found edge state is returned.
Else, a triple `edge, track1, track2` is returned, where `track1` and `track2` are the
tracks along the edge within the basin of attraction of `u1` and `u2`, respectively;
`edge` is the track along the edge derived from `(track1 + track2)/2`.

!!! warning
    May behave erroneously when run with `solver = SimpleATsit5()`, which is the default
    solver for `AttractorsViaProximity`. The recommended solver is `Vern9()`.
"""
function edgetracking(sys::StochSystem, u1::State, u2::State, attractors::Vector;
    eps1=1e-9,
    eps2=1e-8,
    converge=1e-5,
    dt=0.01,
    tmax=Inf,
    ϵ_mapper=0.1,
    solver=Vern9(),
    maxit=100,
    verbose=true,
    output_all=false,
    kwargs...)
    
    diffeq = (;alg = solver)
    cds = to_cds(sys)
    attrs = Dict(i => Dataset([attractors[i]]) for i in 1:length(attractors))
    pinteg = parallel_integrator(cds, [u1, u2]; diffeq)
    mapper = AttractorsViaProximity(cds, attrs, ϵ_mapper; diffeq, dt, kwargs...)
    
    track_edge(pinteg, mapper;
        eps1=eps1, eps2=eps2, converge=converge, dt=dt, tmax=tmax, maxit=maxit,
        verbose=verbose, output_all=output_all)
end;

"""
    bisect_to_edge(sys::StochSystem, u1::State, u2::State, attractors::Vector; kwargs...)
"""
function bisect_to_edge(sys::StochSystem, u1::State, u2::State, attractors::Vector;
    eps1=1e-9,
    dt=0.01,
    ϵ_mapper=0.1,
    solver=Vern9(),
    kwargs...)

    diffeq = (;alg = solver)
    cds = to_cds(sys)
    attrs = Dict(i => Dataset([attractors[i]]) for i in 1:length(attractors))
    pinteg = parallel_integrator(cds, [u1, u2]; diffeq)
    mapper = AttractorsViaProximity(cds, attrs, ϵ_mapper; diffeq, dt, kwargs...)

    bisect_to_edge(pinteg, mapper, abstol=eps1)
end

"""
    track_edge(pinteg, mapper::AttractorMapper; kwargs...)
Internal function for the edge tracking algorithm. See [`edgetracking`](@ref) for details.
"""
function track_edge(pinteg, mapper::AttractorMapper;
    eps1=1e-9, eps2=1e-8, converge=1e-5, dt=0.01, tmax=Inf, maxit=100,
    verbose=false, output_all=false)
    
    verbose && println("=== Starting edge tracking algorithm ===")
    
    u1, u2 = bisect_to_edge(pinteg, mapper; abstol=eps1)
    edgestate = (u1 + u2)/2
    
    if output_all
        track1 = [u1]
        track2 = [u2]
        edge = [edgestate]
    end
    verbose && println("... Iteration 1: Edge at $(edgestate)")

    correction = converge + 1
    counter = 1
    while (correction > converge) & (maxit > counter)
        reinit!(pinteg, [u1,u2])
        state = edgestate
        distance = norm(get_state(pinteg, 1)-get_state(pinteg, 2))
        T = 0
        while (distance < eps2) && (T < tmax)
            step!(pinteg, dt)
            distance = norm(get_state(pinteg, 1)-get_state(pinteg, 2))
            T += dt
        end
        u1, u2 = bisect_to_edge(pinteg, mapper; abstol=eps1)
        edgestate = (u1 + u2)/2
        correction = norm(edgestate - state)
        counter += 1

        if output_all
            push!(track1, u1)
            push!(track2, u2)
            push!(edge, edgestate)
        end
        
        verbose && println("... Iteration $(counter): Edge at $(edgestate)")
        (counter == maxit) && @error("Reached maximum number of $(maxit) iterations; did not converge.")
    end

    (counter < maxit) && println("Edge-tracking converged after $(counter) iterations.")

    if output_all
        return edge, track1, track2
    else
        return edgestate
    end
end;

"""
    bisect_to_edge(pinteg, mapper; abstol=1e-9)

Find the basin boundary between two states `u1, u2 = get_states(pinteg)` by bisecting along a
straight line in phase space. The states `u1` and `u2` must belong to different basins.
Returns two new states located on either side of the basin boundary at a maximum distance of
`abstol` between each other.

`pinteg` is a `parallel_integrator` with two states. The `mapper` must be an `AttractorMapper`
of subtype `AttractorsViaProximity` or `AttractorsViaRecurrences`.

Note: If the straight line between `u1` and `u2` intersects the basin boundary multiple
times, the method will find one of these intersection points. If more than two attractors
exist, one of the two returned states may belong to a different basin than the initial
conditions `u1` and `u2`. A warning is raised if the bisection involves a third basin.

# Keyword arguments
* `abstol = 1e-9`: The maximum (Euclidean) distance between the two returned states.
"""
function bisect_to_edge(pinteg, mapper::AttractorMapper; abstol=1e-9)
    u1, u2 = get_states(pinteg)
    idx1, idx2 = mapper(u1), mapper(u2)
    
    if idx1 == idx2
        error("Both initial conditions belong to the same basin of attraction.")
    end
    
    distance = norm(u1-u2)
    while distance > abstol
        u_new = (u1 + u2)/2
        idx_new = mapper(u_new)
        
        (idx_new == -1) ? error("New point lies exactly on the basin boundary.") : nothing
        
        if idx_new == idx1
            u1 = u_new
        else 
            u2 = u_new
            if idx_new != idx2
                @warn "New bisection point belongs to a third basin of attraction."
            end
        end    
        distance = norm(u1-u2)
    end
    [u1, u2]
end;

"""
    attractor_mapper(sys::StochSystem, attractors, eps=0.01; kwargs...)
Wrapper around the [`AttractorsViaProximity`](https://juliadynamics.github.io/Attractors.jl/dev/attractors/#Attractors.AttractorsViaProximity)
mapper from `DynamicalSystems.jl`. Takes as input `sys` of type [`StochSystem`](@ref), a 
vector `attractors` containing as elements the state vectors of the stable fixed points,
and a parameter `eps` to control the mapping algorithm. For more info, see the
[`docs`](https://juliadynamics.github.io/Attractors.jl/dev/attractors/#Attractors.AttractorsViaProximity).
"""
function attractor_mapper(sys::StochSystem, attractors, eps=0.01;
    dt=0.01, kwargs...)
    AttractorsViaProximity(to_cds(sys), attrs, eps; dt, kwargs...)
end