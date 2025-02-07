module ChaosToolsExt

using CriticalTransitions
using DocStringExtensions
using ChaosTools: ChaosTools, fixedpoints
using DynamicalSystemsBase: CoupledODEs, StateSpaceSet, jacobian
using IntervalArithmetic: IntervalArithmetic, interval

export fixedpoints, intervals_to_box

"""
$(TYPEDSIGNATURES)

Generates a box from specifying the interval limits in each dimension.
* `bmin` (Vector): lower limit of the box in each dimension
* `bmax` (Vector): upper limit

## Example
`intervals_to_box([-2,-1,0], [2,1,1])` returns a 3D box of dimensions `[-2,2] × [-1,1] × [0,1]`.
"""
function CriticalTransitions.intervals_to_box(bmin::Vector, bmax::Vector)
    # Generates a box from specifying the interval limits
    if length(bmin) != length(bmax)
        @warn "bmin and bmax must have the same length."
    end
    return SVector{length(bmin)}([interval(bmin[i], bmax[i]) for i in 1:length(bmin)])
end

"""
$(TYPEDSIGNATURES)

Returns fixed points, their eigenvalues and stability of the system `sys` within the state space volume defined by `bmin` and `bmax`.

> This is a wrapper around the [`fixedpoints`](https://juliadynamics.github.io/ChaosTools.jl/stable/periodicity/#ChaosTools.fixedpoints) function of `DynamicalSystems.jl`.

## Input
* `bmin` (Vector): lower limits of the state space box to be considered, as a vector of coordinates
* `bmax` (Vector): upper limits
* alternatively `box` (IntervalBox) can replace `bmin` and `bmax`

> Example: `fixedpoints(sys, [-2,-1,0], [2,1,1])` finds the fixed points of the 3D system `sys` in a cube defined by the intervals `[-2,2] × [-1,1] × [0,1]`.

## Output
`[fp, eigs, stable]`
* `fp`: `StateSpaceSet` of fixed points
* `eigs`: vector of Jacobian eigenvalues of each fixed point
* `stable`: vector of booleans indicating the stability of each fixed point (`true`=stable, `false`=unstable)

## Additional methods
* `fixedpoints(sys::CoupedSDEs, box)`
"""
function ChaosTools.fixedpoints(sys::CoupledSDEs, bmin::Vector, bmax::Vector)
    box = CriticalTransitions.intervals_to_box(bmin, bmax)
    return fixedpoints(sys::CoupledSDEs, box)
end

function ChaosTools.fixedpoints(sys::CoupledSDEs, box)
    ds = CoupledODEs(sys)
    return fixedpoints(CoupledODEs(sys), box, jacobian(ds))
end

end
