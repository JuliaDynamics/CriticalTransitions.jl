module AttractorsExt

using CriticalTransitions
using Attractors:
    Attractors, AttractorsViaProximity, AttractorMapper
import Attractors: edgetracking, bisect_to_edge
using DynamicalSystemsBase: ParallelDynamicalSystem

export edgetracking, bisect_to_edge

"""
    edgetracking(ds::CoupledSDEs, attractors::Dict; diffeq, kwars...)

Runs the edge tracking algorithm for the deterministic part of the CoupledSDEs `ds`.

## Keyword arguments
- `diffeq=(;alg = Vern9(), reltol=1e-11)`: ODE solver settings
- `kwargs...`: all keyword arguments of `Attractors.edgetracking`
"""
function edgetracking(ds::CoupledSDEs, attractors::Dict;
    diffeq=(;alg = Vern9(), reltol=1e-11),
    kwargs...)
    return Attractors.edgetracking(CoupledODEs(ds; diffeq), attractors; kwargs...)
end

"""
    bisect_to_edge(ds::CoupledSDEs, attractors::Dict;
        u1, u2, bisect_thresh, diffeq, verbose, kwargs...)

Runs the `bisect_to_edge` function for the deterministic part of the CoupledSDEs `ds`.

## Keyword arguments
- `u1`: State 1 (default: first state in `attractors`)
- `u2`: State 2 (default: second state in `attractors`)
- `bisect_thresh=1e-6`: distance threshold
- `diffeq=(;alg = Vern9(), reltol=1e-11)`: ODE solver settings
- `verbose=false`: Verbosity of output
- `ϵ_mapper=nothing`: ϵ argument of `Attractors.AttractorsViaProximity`
- `kwargs...`: Keyword arguments passed to `Attractors.AttractorsViaProximity`
"""
function bisect_to_edge(ds::CoupledSDEs, attractors::Dict;
    u1=collect(values(attractors))[1][1],
    u2=collect(values(attractors))[2][1],
    bisect_thresh=1e-6,
    diffeq=(;alg = Vern9(), reltol=1e-11),
    ϵ_mapper=nothing,
    verbose=false,
    kwargs...)

    odes = CoupledODEs(ds; diffeq)
    pds = Attractors.ParallelDynamicalSystem(odes, [u1, u2])
    mapper = Attractors.AttractorsViaProximity(odes, attractors, ϵ_mapper; kwargs...)

    return Attractors.bisect_to_edge(pds, mapper; bisect_thresh, verbose)
end

end # module AttractorsExt
