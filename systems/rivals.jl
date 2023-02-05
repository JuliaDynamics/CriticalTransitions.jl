"""
    rivals!(du, u, p, t)
In-place definition of the Rivals system.

See also [`rivals`](@ref).
"""
function rivals!(du, u, p, t)
    x, y = u
    ϵ, α₁, α₂, β₁, β₂ = p[1]

    du[1] = x(x-α₁)(1-x)-β₁x*y
    du[2] = ϵ*(y*(y-α₂)*(1-y)-β₂x*y)

end

"""
    rivals(u, p, t)
Out-of-place definition of the Rivals system.

See also [`rivals!`](@ref).
"""
function rivals(u, p, t)
    x, y = u    
    ϵ, α₁, α₂, β₁, β₂ = p[1]

    dx = x(x-α₁)(1-x)-β₁x*y
    dy = ϵ*(y*(y-α₂)*(1-y)-β₂x*y)

    SVector{2}(dx, dy)

end