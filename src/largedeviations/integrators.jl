# abstract type AbstractIntegrator end

# struct Euler{T} <: AbstractIntegrator
#     ϵ::T
#     integ
# end

# struct RK4{T} <: AbstractIntegrator
#     ϵ::T
# end

"""
$(TYPEDSIGNATURES)

Perform an in-place Euler step on the path `x` with the gradient `sys.H_p` of the extented phase space Hamiltonian with step size `ϵ`.

"""
function updateEuler!(x::Matrix, sys::SgmamSystem, ϵ::Real)
    return x += ϵ * sys.H_p(x, 0 * x) # euler integration
end
"""
$(TYPEDSIGNATURES)

Perform an in-place Euler step on the path `x` with the gradient `sys.H_p` of the extented phase space Hamiltonian with step size `ϵ`.

"""
function updateEuler!(x::StateSpaceSet, sys::SgmamSystem, ϵ::Real)
    return x += ϵ .* vec(sys.H_p(x, 0 * Matrix(x))) # euler integration
end
"""
$(TYPEDSIGNATURES)

Perform an in-place Euler step on the path `x` with the gradient `b` of the extented phase space Hamiltonian with step size `ϵ`.

"""
function updateEuler!(x::StateSpaceSet, b::Function, ϵ::Real)
    # do not touch the initial and final points to allow for string computation
    # between points that are not stable fixed points
    for j in 2:(length(x) - 1)
        x[j] += ϵ * b(x[j]) # euler integration
    end
end
"""
$(TYPEDSIGNATURES)

Perform an in-place Euler step on the path `x` with the gradient `b` of the extented phase space Hamiltonian with step size `ϵ`.

"""
function updateEuler!(x::Matrix, b::Function, ϵ::Real)
    # do not touch the initial and final points to allow for string computation
    # between points that are not stable fixed points
    for j in 2:(size(x)[2] - 1)
        x[:, j] += ϵ * b(x[:, j]) # euler integration
    end
end

"""
$(TYPEDSIGNATURES)

Perform an in-place Runga-Kuta step of order 4 on the path `x` with the gradient `sys.H_p` of the extented phase space Hamiltonian with step size `ϵ`.

"""

function updateRK4!(x::Matrix, sys::SgmamSystem, ϵ::Real)
    p = 0 * x
    gradV1 = sys.H_p(x, p)
    gradV2 = sys.H_p(x - 0.5 * ϵ * gradV1, p)
    gradV3 = sys.H_p(x - 0.5 * ϵ * gradV2, p)
    gradV4 = sys.H_p(x - ϵ * gradV3, p)
    return x += ϵ / 6 * (gradV1 + 2 * gradV2 + 2 * gradV3 + gradV4) # ±?
end
"""
$(TYPEDSIGNATURES)

Perform an in-place Euler step on the path `x` with the gradient `sys.H_p` of the extented phase space Hamiltonian with step size `ϵ`.

"""
function updateRK4!(x::StateSpaceSet, sys::SgmamSystem, ϵ::Real)
    p = 0 * Matrix(x)
    gradV1 = vec(sys.H_p(x, p))
    gradV2 = vec(sys.H_p(x - 0.5 * ϵ * gradV1, p))
    gradV3 = vec(sys.H_p(x - 0.5 * ϵ * gradV2, p))
    gradV4 = vec(sys.H_p(x - ϵ * gradV3, p))
    return x += ϵ / 6 * (gradV1 + 2 * gradV2 + 2 * gradV3 + gradV4) # ±?
end
"""
$(TYPEDSIGNATURES)

Perform an in-place Euler step on the path `x` with the gradient `b` of the extented phase space Hamiltonian with step size `ϵ`.

"""
function updateRK4!(x::StateSpaceSet, b::Function, ϵ::Real)
    # do not touch the initial and final points to allow for string computation
    # between points that are not stable fixed points
    for j in 2:(length(x) - 1)
        gradV1 = b(x[j])
        gradV2 = b(x[j] - 0.5 * ϵ * gradV1)
        gradV3 = b(x[j] - 0.5 * ϵ * gradV2)
        gradV4 = b(x[j] - ϵ * gradV3)
        x[j] += ϵ / 6 * (gradV1 + 2 * gradV2 + 2 * gradV3 + gradV4) # ±?
    end
end
"""
$(TYPEDSIGNATURES)

Perform an in-place Euler step on the path `x` with the gradient `b` of the extented phase space Hamiltonian with step size `ϵ`.

"""
function updateRK4!(x::Matrix, b::Function, ϵ::Real)
    # do not touch the initial and final points to allow for string computation
    # between points that are not stable fixed points
    for j in 2:(size(x)[2] - 1)
        gradV1 = b(x[:, j])
        gradV2 = b(x[:, j] - 0.5 * ϵ * gradV1)
        gradV3 = b(x[:, j] - 0.5 * ϵ * gradV2)
        gradV4 = b(x[:, j] - ϵ * gradV3)
        x[:, j] += ϵ / 6 * (gradV1 + 2 * gradV2 + 2 * gradV3 + gradV4) # ±?
    end
end
