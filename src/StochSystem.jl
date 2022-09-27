using LinearAlgebra

include("utils.jl")

# Define custom types
Parameters = Union{Vector{Float64}, Vector{Int64}, Nothing};
CovMatrix = Union{Matrix, UniformScaling{Bool}, Diagonal{Bool, Vector{Bool}}};

# Define StochSystem
struct StochSystem
    f::Function
    pf::Parameters
    σ::Float64
    g::Function
    pg::Parameters
    Σ::CovMatrix
    process::Any
end

# Additional methods of StochSystem
StochSystem(f, pf) = StochSystem(f, pf, 1.0, idfunc, nothing, I, "white-gauss")
StochSystem(f, pf, σ) = StochSystem(f, pf, σ, idfunc, nothing, I, "white-gauss")
StochSystem(f, pf, σ, Σ) = StochSystem(f, pf, σ, idfunc, nothing, Σ, "white-gauss")

function σg(sys::StochSystem)
    # Multiplies the noise strength σ of a StochSystem with its noise function g
    if is_iip(sys.f)
        g_iip(du,u,p,t) = sys.σ .* sys.g(du,u,p,t)
    else
        g_oop(u,p,t) = sys.σ .* sys.g(u,p,t)
    end
end;

function p(sys::StochSystem)
    # Concatenates the deterministic and stochastic parameter vectors pf and pg of a StochSystem
    [sys.pf, sys.pg]
end;
