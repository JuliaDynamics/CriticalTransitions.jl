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
function anorm(vec, A; square = false)
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
function subnorm(vec; directions = 1:length(vec))
    sum = 0
    for i in directions
        sum += vec[i]^2
    end
    return sqrt(sum)
end;

"""
$(TYPEDSIGNATURES)

Returns ``\\mathbf{Q}\\cdot D/\\mathrm{tr}(\\mathbf{Q})``, the trace-normalized covariance
(``\\mathrm{tr} = D``, average eigenvalue 1).

For an SDE built as ``\\mathrm{d}\\mathbf{x} = \\mathbf{b}\\,\\mathrm{d}t + \\sigma\\sqrt{\\mathbf{Q}}\\,\\mathrm{d}\\mathbf{W}``
the SDE only fixes the product ``\\sigma^2\\mathbf{Q} = `` `covariance_matrix(sys)`; this
function picks the canonical pair ``(\\sigma_{\\mathrm{eff}}, \\mathbf{Q}_{\\mathrm{can}})``
defined by ``\\mathrm{tr}(\\mathbf{Q}_{\\mathrm{can}}) = D`` (see [`noise_strength`](@ref) for the
matching ``\\sigma_{\\mathrm{eff}}``). Used internally by [`fw_action`](@ref),
[`om_action`](@ref), and [`geometric_action`](@ref). For diagonal positive ``\\mathbf{Q}``
this coincides numerically with dividing by ``L_1(\\mathbf{Q})/D``.

See [Large deviation theory](@ref) for the rotation-invariance and conversion factors.
"""
function normalize_covariance!(covariance)
    D = size(covariance, 1)
    return covariance * D / tr(covariance)
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
