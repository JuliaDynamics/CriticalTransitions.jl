export edgetracking
using DynamicalSystems, DynamicalSystemsBase

function norm(x)
    sum = 0
    for i in 1:length(x)
        sum += x[i]^2
    end
    sqrt(sum)
end

function edgetracking(sys::StochSystem, u1, u2, attractors::Dict;
    dist_bisect=1e-9,
    diverge_factor=10,
    converge=1e-5,
    Δt=1)
    
    ds = to_cds(sys)
    pinteg = parallel_integrator(ds, [u1, u2])
    mapper = AttractorsViaProximity(ds, attractors)
    
    track_edge(pinteg, mapper;
        dist_bisect=dist_bisect, diverge_factor=diverge_factor, converge=converge, Δt=Δt)
end

function track_edge(sys::StochSystem, mapper, init1, init2; dist_bisect=1e-9, diverge_factor=2, converge=1e-5, Δt=0.01)
    println("+++++++ Start edge tracking algorithm ++++++++++")
    u1, u2 = bisect_to_edge(mapper, init1, init2; abstol=dist_bisect)
    edgestate = (u1 + u2)/2
    println("u1", mapper(u1))
    println("u2", mapper(u2))
    
    correction = converge + 1
    counter = 0
    while correction > converge
        state = edgestate
        println("counter ", counter)
        println(mapper(u1), u1)
        println(mapper(u2), u2)
        #reinit!(pinteg, [u1, u2])
        distance = norm(u1-u2)
        ii = 0
        while distance < dist_bisect*diverge_factor
            x1, x2 = u1, u2
            println("ii ", ii, get_states(pinteg))
            ii += 1
            x1 = relax(sys, x1; dt=Δt, tmax=Δt)
            step!(pinteg, Δt)
            distance = norm(get_state(pinteg, 1)-get_state(pinteg, 2))
        end
        println(get_states(pinteg))
        println(mapper(get_state(pinteg, 1)), mapper(get_state(pinteg, 2)))
        u1, u2 = bisect_to_edge(pinteg, mapper; abstol=dist_bisect)
        edgestate = (u1 + u2)/2
        
        correction = norm(edgestate - state)
        counter += 1
    end
    edgestate
end

"""
    bisect_to_edge(mapper, u1, u2; abstol=1e-9)

Find the basin boundary between two states `u1`, `u2` by bisecting along a
straight line in phase space. The states `u1` and `u2` must belong to different basins.
Returns two new states located on either side of the basin boundary at a distance of less
than `abstol` between each other.

The `mapper` must be an `AttractorMapper` of subtype `AttractorsViaProximity` or
`AttractorsViaRecurrences`.

Note: If the straight line between `u1` and `u2` intersects the basin boundary multiple
times, the method will find one of these intersection points. If more than two attractors
exist, one of the two returned states may belong to a different basin than the initial
conditions `u1` and `u2`. A warning is raised if the bisection involves a third basin.

# Keyword arguments
* `abstol = 1e-9`: The algorithm stops once the two updated states are less than `abstol`
apart from each other, measured in Euclidean distance in state space.
"""
function bisect_to_edge(mapper::AttractorMapper, u1, u2; abstol=1e-9)
    #u1, u2 = get_states(pinteg)
    idx1, idx2 = mapper(u1), mapper(u2)
    
    if idx1 == idx2
        error("Both initial conditions belong to the same basin of attraction.")
    end
    
    distance = norm(u1-u2)
    while distance > abstol
        u_new = (u1 + u2)/2
        idx_new = mapper(u_new)
        # todo: what if u_new lies on boundary
        
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
end


using CriticalTransitions

p = [1.,3.,1.,1.,1.,0.]
sys = StochSystem(FitzHughNagumo, p, 2)

u1 = [0.2801377011835575, 1.0]
u2 = [0.27013770025223494, 1.0]

sim1 = relax(sys, u1; dt=0.01, tmax=100)
sim2 = relax(sys, u2; dt=0.01, tmax=100)

sim1[:,end]
sim2[:,end]