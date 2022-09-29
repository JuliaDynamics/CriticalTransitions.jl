using StaticArrays

function FitzHughNagumo!(du, u, p, t)
    x, y = u
    ϵ, β, α, γ, κ, I = p[1]

    du[1] = (-α*x^3 + γ*x - κ*y + I)/ϵ
    du[2] = -β*y + x
end

function FitzHughNagumo(u,p,t)
    x, y = u
    ϵ, β, α, γ, κ, I = p[1]

    dx = (-α*x^3 + γ*x - κ*y + I)/ϵ
    dy = -β*y + x

    SVector{2}([dx, dy])
end