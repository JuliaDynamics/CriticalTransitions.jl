# Define custom types
Parameters = Union{Vector{Any}, Nothing};
CovMatrix = Union{Matrix, UniformScaling{Bool}, Diagonal{Bool, Vector{Bool}}};
State = Union{Vector, SVector}

# Define StochSystem
"""
    StochSystem(f, pf, dim, σ, g, pg, Σ, process)
Defines a stochastic dynamical system in `CriticalTransitions.jl`. See [documentation](https://reykboerner.github.io/CriticalTransitions.jl/dev/man/stochsystem/).
"""
struct StochSystem
    f::Function
    pf::Parameters
    u::State
    σ::Float64
    g::Function
    pg::Parameters
    Σ::CovMatrix
    process::Any
end;

# Methods of StochSystem
StochSystem(f, pf, u) = StochSystem(f, pf, u, 0.0, idfunc, nothing, I(length(u)), "WhiteGauss")
StochSystem(f, pf, u, σ) = StochSystem(f, pf, u, σ, idfunc, nothing, I(length(u)), "WhiteGauss")
StochSystem(f, pf, u, σ, Σ) = StochSystem(f, pf, u, σ, idfunc, nothing, Σ, "WhiteGauss")

# Core functions acting on StochSystem
"""
    drift(sys::StochSystem, x::State)
Returns the drift field ``b(x)`` of the StochSystem `sys` at the state vector `x`.
"""
drift(sys::StochSystem, x::State) = sys.f(x, [sys.pf], 0)

"""
    σg(sys::StochSystem)
Multiplies the noise strength `σ` of a StochSystem `sys` with its noise function `g`.
"""
function σg(sys::StochSystem)
    if is_iip(sys.f)
        g_iip(du,u,p,t) = sys.σ .* sys.g(du,u,p,t)
    else
        g_oop(u,p,t) = sys.σ .* sys.g(u,p,t)
    end
end

"""
    p(sys::StochSystem)
Concatenates the deterministic and stochastic parameter vectors `pf` and `pg` of a StochSystem `sys`.
"""
function p(sys::StochSystem)
    [sys.pf, sys.pg]
end;

"""
    CoupledODEs(sys::StochSystem; diffeq, t0=0.0)
Converts a [`StochSystem`](@ref) into [`CoupledODEs`](https://juliadynamics.github.io/DynamicalSystems.jl/stable/tutorial/#DynamicalSystemsBase.CoupledODEs)
from DynamicalSystems.jl.
"""
CoupledODEs(sys::StochSystem; diffeq=DynamicalSystemsBase.DEFAULT_DIFFEQ, t0=0.0) =
DynamicalSystems.CoupledODEs(sys.f, SVector{2}(sys.u), [sys.pf]; diffeq=diffeq, t0=t0)

to_cds(sys::StochSystem) = CoupledODEs(sys)

"""
    StochSystem(ds::CoupledODEs, σ, g, pg, Σ, process)
Converts a [`CoupledODEs`](https://juliadynamics.github.io/DynamicalSystems.jl/stable/tutorial/#DynamicalSystemsBase.CoupledODEs)
system into a [`StochSystem`](@ref).
"""
StochSystem(ds::DynamicalSystemsBase.CoupledODEs, σ=0.0, g=idfunc, pg=nothing, Σ=I(length(get_state(ds))), process="WhiteGauss") =
StochSystem(dynamic_rule(ds), [ds.p0], get_state(ds), σ, g, pg, Σ, process)