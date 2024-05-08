"""
Utility functions for CriticalTransitions.jl
Functions in this file are independent of types/functions defined in CriticalTransitions
"""

"""
$(TYPEDSIGNATURES)

Identity function for noise function `StochSystem.g` (out-of-place).
"""
function idfunc(u, p, t)
    ones(length(u))
end;

"""
$(TYPEDSIGNATURES)

Identity function for noise function `StochSystem.g` (in-place).
"""
function idfunc!(du, u, p, t)
    du .= ones(length(u))
    nothing
end;

"""
$(TYPEDSIGNATURES)

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
    du[idx] = 1.
    nothing
end;

function additive_idx(u,p,t,idx)
    du = zeros(length(u))
    du[idx] = 1.
    SVector{length(u)}(du)
end;

function multiplicative_idx!(du, u, p, t, idx)
    du[idx] = u[idx]
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
$(TYPEDSIGNATURES)

Calculates the generalized ``A``-norm of the vector `vec`,
``||v||_A := \\sqrt(v^\\top \\cdot A \\cdot v)``,
where `A` is a square matrix of dimension `(lenth(vec) x length(vec))`.

# Keyword arguments
* `square`: if `true`, returns square of norm; else, returns norm.
"""
function anorm(vec, A; square=false)
    normsquared = dot(vec, A, vec)
    square ? normsquared : sqrt(normsquared)
end;

"""
$(TYPEDSIGNATURES)

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

# Central finite difference, second derivative
function central(f, idx, dz)
    (f[idx+1] - f[idx-1])/(2*dz)
end;

function central2(f, idx, dz)
    (f[idx+1] - 2f[idx] + f[idx-1])/(dz^2)
end;

"""
$(TYPEDSIGNATURES)
Smooth approximation of `abs(x)`, ``|x| = x \\tanh(\\xi x)``, where ``xi`` controls the
accuracy of the approximation. The exact absolute value function is obtained in the limit
``\\xi \\to \\infty``.
"""
function smoothabs(x, xi=1000)
    x*tanh(x*xi)
end;

"""
"""
function diag_noise_funtion(σ; in_place = false)
    if in_place
        return (du, u, p, t) -> σ .* idfunc!(du, u, p, t)
    else
        return (u, p, t) -> σ .* idfunc(u, p, t)
    end
end
function diag_noise_funtion(σ, g)
    if SciMLBase.isinplace(g, 4)
        return (du, u, p, t) -> σ .* g(du, u, p, t)
    else
        return (u, p, t) -> σ .* g(u, p, t)
    end
end
