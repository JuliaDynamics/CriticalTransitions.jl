"""
Utility functions for CriticalTransitions.jl
Functions in this file are independent of types/functions defined in CriticalTransitions
"""

"""
    idfunc(u, p, t)
Identity function for noise function `StochSystem.g` (out-of-place).
"""
function idfunc(u, p, t)
    ones(length(u))
end;

"""
    idfunc!(du, u, p, t)
Identity function for noise function `StochSystem.g` (in-place).
"""
function idfunc!(du, u, p, t)
    du = ones(length(u))
end;

"""
    is_iip(f::Function)
Asserts if f is in-place (true) or out-of-place (false).

> Warning: This function simply checks if there is a `!` in the function name. Thus, if you do not add a `!` to in-place function names (as recommended by the Julia style guide), this test will not work.
"""
function is_iip(f::Function)
    occursin("!", String(Symbol(f)))
end;

"""
    intervals_to_box(bmin::Vector, bmax::Vector)
Generates a box from specifying the interval limits in each dimension.
* `bmin` (Vector): lower limit of the box in each dimension
* `bmax` (Vector): upper limit

## Example 
`intervals_to_box([-2,-1,0], [2,1,1])` returns a 3D box of dimensions `[-2,2] × [-1,1] × [0,1]`.
"""
function intervals_to_box(bmin::Vector, bmax::Vector)
    # Generates a box from specifying the interval limits
    intervals = []
    dim = length(bmin)
    for i in 1:dim
        push!(intervals, bmin[i]..bmax[i])
    end
    box = intervals[1]
    for i in 2:dim
        box = box × intervals[i]
    end
    box
end;