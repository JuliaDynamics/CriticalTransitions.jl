module ChaosToolsExt

using CriticalTransitions
using DocStringExtensions
using ForwardDiff
using ChaosTools: ChaosTools, fixedpoints

export fixedpoints

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
* `fp`: `Dataset` of fixed points
* `eigs`: vector of Jacobian eigenvalues of each fixed point
* `stable`: vector of booleans indicating the stability of each fixed point (`true`=stable, `false`=unstable)

## Additional methods
* `fixedpoints(sys::CoupedSDEs, box)`
"""
function ChaosTools.fixedpoints(sys::CoupledSDEs, bmin::Vector, bmax::Vector)
    box = intervals_to_box(bmin, bmax)
    return fixedpoints(sys::CoupledSDEs, box)
end

function ChaosTools.fixedpoints(sys::CoupledSDEs, box)
    jac(u, p, t) = ForwardDiff.jacobian((x) -> sys.integ.f(x, p, t), u)
    return fixedpoints(CoupledODEs(sys), box, jac)
end

# function saddles_idx(fps::Tuple)
#     num = size(fps[1],1); # number of fixed points
#     dim = size(fps[1],2); # dimension of the system
#     eigenvalues = fps[2];
#     idx = [false for i ∈ 1:num];
#     for ii ∈ 1:num
#         imag_parts = [imag(eigenvalues[ii][jj]) for jj ∈ 1:dim]
#         if all(imag_parts.==0) # we have purely real eigenvalues
#             real_parts = [real(eigenvalues[ii][jj]) for jj ∈ 1:dim];
#             if prod(real_parts) < 0 # we have at least positive eigenvalue and at least one negative eigenvalue
#                 idx[ii] = true;
#             end
#         end
#     end
#     idx
# end

# function repellers_idx(fps::Tuple)
#     num = size(fps[1],1); # number of fixed points
#     dim = size(fps[1],2); # dimension of the system
#     eigenvalues = fps[2];
#     idx = [false for i ∈ 1:num];
#     for ii ∈ 1:num
#         real_parts = [real(eigenvalues[ii][jj]) for jj ∈ 1:dim];
#         if all(real_parts .> 0)
#             idx[ii] = true;
#         end
#     end
#     idx
# end

# function attractors_idx(fps::Tuple)
#     num = size(fps[1],1); # number of fixed points
#     dim = size(fps[1],2); # dimension of the system
#     eigenvalues = fps[2];
#     idx = [false for i ∈ 1:num];
#     for ii ∈ 1:num
#         real_parts = [real(eigenvalues[ii][jj]) for jj ∈ 1:dim];
#         if all(real_parts .< 0)
#             idx[ii] = true;
#         end
#     end
#     idx
# end

end
