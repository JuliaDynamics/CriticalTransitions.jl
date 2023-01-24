"""
Dynamical systems specification file
"""

# FitzHughNagumo model

"""
    FitzHughNagumo!(du, u, p, t)
In-place definition of the FitzHugh-Nagumo system.

See also [`FitzHughNagumo`](@ref).
"""
function FitzHughNagumo!(du, u, p, t)
    x, y = u
    ϵ, β, α, γ, κ, I = p[1]

    du[1] = (-α*x^3 + γ*x - κ*y + I)/ϵ
    du[2] = -β*y + x
end

"""
    FitzHughNagumo(u, p, t)
Out-of-place definition of the FitzHugh-Nagumo system.

See also [`FitzHughNagumo!`](@ref).
"""
function FitzHughNagumo(u,p,t)
    x, y = u
    ϵ, β, α, γ, κ, I = p[1]

    dx = (-α*x^3 + γ*x - κ*y + I)/ϵ
    dy = -β*y + x

    SVector{2}([dx, dy])
end