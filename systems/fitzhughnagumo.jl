"""
Dynamical systems specification file
"""

"""
    fitzhugh_nagumo!(du, u, p, t)
In-place definition of the FitzHugh-Nagumo system.

See also [`fitzhugh_nagumo`](@ref).
"""
function fitzhugh_nagumo!(du, u, p, t)
    x, y = u
    ϵ, β, α, γ, κ, I = p

    du[1] = (-α * x^3 + γ * x - κ * y + I) / ϵ
    du[2] = -β * y + x
    return nothing
end

"""
    fitzhugh_nagumo(u, p, t)
Out-of-place definition of the FitzHugh-Nagumo system.

See also [`fitzhugh_nagumo!`](@ref).
"""
function fitzhugh_nagumo(u, p, t)
    x, y = u
    ϵ, β, α, γ, κ, I = p

    dx = (-α * x^3 + γ * x - κ * y + I) / ϵ
    dy = -β * y + x

    return SA[dx, dy]
end