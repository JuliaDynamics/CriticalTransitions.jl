function FitzHughNagumo!(du, u, p, t)
    x, y = u
    ϵ, β, α, γ, κ, I = p

    du[1] = (-α*x^3 + γ*x - κ*y + I)/ϵ
    du[2] = -β*y + x
end

function FitzHughNagumo(u,p,t)
    x, y = u
    ϵ, β, α, γ, κ, I = p

    dx = (-α*x^3 + γ*x - κ*y + I)/ϵ
    dy = -β*y + x

    [dx, dy]
end