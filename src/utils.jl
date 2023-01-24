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

function additive_idx!(du, u, p, t, idx)
    du[indx] .= 1.
end;

function additive_idx(u,p,t,idx)
    du = zeros(length(u))
    du[idx] .= 1.
    SVector{length(u)}(du)
end;

function multiplicative_idx!(du, u, p, t, idx)
    du[indx] .= u[idx]
end;

function multiplicative_idx(u, p, t, idx)
    du = zeros(length(u))
    du[idx] = u[idx]
    SVector{length(u)}(du)
end


function multiplicative_first!(du, u, p, t)
    du[1] = u[1];
end;

function multiplicative_first(u, p, t)
    
    du = zeros(length(u));
    du[1] = u[1];

    SVector{length(u)}(du)
end;

function additive_first!(du, u, p, t)
    du[1] = 1;
end;

function additive_first(u, p, t)
    
    du = zeros(length(u));
    du[1] = 1;

    SVector{length(u)}(du)
end;

"""
    anorm(vec, A; square=false)
Calculates the norm of the vector `vec` with respect to the A-metric, where `A` is a 
square matrix of dimension `(lenth(vec) x length(vec))`.

# Keyword arguments
* `square`: if `true`, returns square of norm; else, returns norm.
"""
function anorm(vec, A; square=false)
    normsquared = dot(vec, A * vec)
    if square
        return normsquared
    else
        return sqrt(normsquared)
    end
end;

"""
    subnorm(vec; directions=[1,...,N])
Returns the Euclidean norm of the vector `vec`; however, if `directions` are specified, the
norm is only computed along the specified components of the vector, i.e. in a subspace of
dimension `length(directions)`.

# Keyword arguments
* `directions`: a vector containing the indices of the components of `vec` to
be included. Defaults to `1:length(vec)`, i.e. all components.

# Example
`subnorm([3,7,4]; directions=[1,3]` calculates the norm of only the 1st and 3rd components
of the vector [3,7,4]:

``\\sqrt{3^2+4^2} = 5``.
"""
function subnorm(vec; directions=1:length(vec))
    sum = 0
    for i in directions
        sum += vec[i]^2
    end
    sqrt(sum)
end;