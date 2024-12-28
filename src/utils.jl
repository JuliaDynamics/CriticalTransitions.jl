"""
Utility functions for CriticalTransitions.jl
Functions in this file are independent of types/functions defined in CriticalTransitions
"""

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

"""
$(TYPEDSIGNATURES)

Normalizes the covariance matrix ``Q`` (in-place) by dividing it by
``L_1(Q)/\\text{dim}(Q)``, where ``L_1(Q)`` is the L1 matrix norm of ``Q`` and
``\\text{dim}(Q)`` is the dimension of ``Q``.
"""
function normalize_covariance!(covariance)
    l1norm = norm(covariance, 1)
    dim = size(covariance)[1]
    return covariance * dim / l1norm
end

# Central finite difference, second derivative
function central(f, idx, dz)
    return (f[idx + 1] - f[idx - 1]) / (2 * dz)
end;

function central2(f, idx, dz)
    return (f[idx + 1] - 2f[idx] + f[idx - 1]) / (dz^2)
end;

function Base.diff(a::StateSpaceSet)
    r = length(a)
    r0 = 1:(r - 1)
    r1 = 2:r

    return StateSpaceSet([a[r1[i]] - a[r0[i]] for i in 1:length(r0)])
end
