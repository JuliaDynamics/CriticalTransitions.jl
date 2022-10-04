using LinearAlgebra, StaticArrays

include("utils.jl")

# Define custom types
Parameters = Union{Vector{Float64}, Vector{Int64}, Nothing};
CovMatrix = Union{Matrix, UniformScaling{Bool}, Diagonal{Bool, Vector{Bool}}};
State = Union{Vector, SVector}

# Define StochSystem
struct StochSystem
    f::Function
    pf::Parameters
    dim::Int64
    σ::Float64
    g::Function
    pg::Parameters
    Σ::CovMatrix
    process::Any
end;

# Methods of StochSystem
StochSystem(f, pf, dim) = StochSystem(f, pf, dim, 1.0, idfunc, nothing, I(dim), "WhiteGauss")
StochSystem(f, pf, dim, σ) = StochSystem(f, pf, dim, σ, idfunc, nothing, I(dim), "WhiteGauss")
StochSystem(f, pf, dim, σ, Σ) = StochSystem(f, pf, dim, σ, idfunc, nothing, Σ, "WhiteGauss")

# Core functions acting on StochSystem
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

function tocds(sys::StochSystem; state=zeros(sys.dim))
    # Turns a StochSystem into a ContinuousDynamicalSystem
    ContinuousDynamicalSystem(sys.f, state, p(sys))
end;