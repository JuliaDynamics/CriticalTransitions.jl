using LinearAlgebra

include("utils.jl")

Parameters = Union{Vector{Float64}, Vector{Int64}, Nothing};
CovMatrix = Union{Matrix, UniformScaling{Bool}, Diagonal{Bool, Vector{Bool}}};

struct StochSystem
    f::Function
    p::Parameters
    σ::Float64
    g::Function
    pg::Parameters
    Σ::CovMatrix
    process::Any
end

StochSystem(f, p) = StochSystem(f, p, 1.0, idfunc, nothing, I, "white-gauss")
StochSystem(f, p, σ) = StochSystem(f, p, σ, idfunc, nothing, I, "white-gauss")