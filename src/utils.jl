"""
Utility functions for CriticalTransitions.jl
Functions in this file are independent of types/functions defined in CriticalTransitions
"""

"""
$(TYPEDSIGNATURES)

Identity function for a diffusion function `g` of `CoupledSDEs` (out-of-place).
Equivalent to `(u, p, t) -> ones(length(u))`,
"""
function idfunc(u, p, t)
    return typeof(u)(ones(eltype(u), length(u)))
end;

"""
$(TYPEDSIGNATURES)

Identity function for a diffusion function `g` of `CoupledSDEs` (in-place).
Equivalent to `idfunc!(du, u, p, t) = (du .= ones(length(u)); return nothing)`
"""
function idfunc!(du, u, p, t)
    du .= ones(eltype(u), length(u))
    return nothing
end;

function σg(σ, g)
    return (u, p, t) -> σ .* g(u, p, t)
end

function σg!(σ, g!)
    function (du, u, p, t)
        g!(du, u, p, t)
        du .*= σ
        return nothing
    end
end

function add_noise_strength(σ, g, IIP)
    newg = IIP ? σg!(σ, g) : σg(σ, g)
    return newg
end

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
        push!(intervals, interval(bmin[i], bmax[i]))
    end
    box = intervals[1]
    for i in 2:dim
        box = IntervalArithmetic.cross(box, intervals[i])
    end
    return box
end;

"""
Calculates the generalized ``A``-norm of the vector `vec`,
``||v||_A := \\sqrt(v^\\top \\cdot A \\cdot v)``,
where `A` is a square matrix of dimension `(length(vec) x length(vec))`.

# Keyword arguments
* `square`: if `true`, returns square of norm; else, returns norm.
"""
function anorm(vec, A; square=false)
    normsquared = dot(vec, A, vec)
    return square ? normsquared : sqrt(normsquared)
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
`subnorm([3,7,4]; directions=[1,3])` calculates the norm of only the 1st and 3rd components
of the vector [3,7,4]:

``\\sqrt{3^2+4^2} = 5``.
"""
function subnorm(vec; directions=1:length(vec))
    sum = 0
    for i in directions
        sum += vec[i]^2
    end
    return sqrt(sum)
end;

# Central finite difference, second derivative
function central(f, idx, dz)
    return (f[idx + 1] - f[idx - 1]) / (2 * dz)
end;

function central2(f, idx, dz)
    return (f[idx + 1] - 2f[idx] + f[idx - 1]) / (dz^2)
end;

"""
$(TYPEDSIGNATURES)
Smooth approximation of `abs(x)`, ``|x| = x \\tanh(\\xi x)``, where ``xi`` controls the
accuracy of the approximation. The exact absolute value function is obtained in the limit
``\\xi \\to \\infty``.
"""
function smoothabs(x, xi=1000)
    return x * tanh(x * xi)
end;
