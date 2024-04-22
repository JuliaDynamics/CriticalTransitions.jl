using DynamicalSystemsBase: CoupledODEs, isinplace, __init, SciMLBase, correct_state
using StochasticDiffEq: SDEProblem, SDEIntegrator
using StochasticDiffEq: EM, SDEProblem

const DEFAULT_SOLVER = SOSRI()
const DEFAULT_DIFFEQ_KWARGS = (abstol = 1e-6, reltol = 1e-6)
const DEFAULT_DIFFEQ = (alg = DEFAULT_SOLVER, DEFAULT_DIFFEQ_KWARGS...)

# Function from user `@xlxs4`, see
# https://github.com/JuliaDynamics/DynamicalSystemsBase.jl/pull/153
_delete(a::NamedTuple, s::Symbol) = NamedTuple{filter(≠(s), keys(a))}(a)
function _decompose_into_solver_and_remaining(diffeq)
    if haskey(diffeq, :alg)
        return (diffeq[:alg], _delete(diffeq, :alg))
    else
        return (DEFAULT_SOLVER, diffeq)
    end
end

# Define custom types
# Parameters = Union{Vector{Any}, Nothing};
# CovMatrix = Union{Matrix, UniformScaling{Bool}, Diagonal{Bool, Vector{Bool}}};
# State = Union{Vector, SVector}

# Define StochSystem
"""
    StochSystem(f, pf, dim, σ, g, pg, Σ, process)

Defines a stochastic dynamical system in `CriticalTransitions.jl`. See [documentation](https://juliadynamics.github.io/CriticalTransitions.jl/dev/man/stochsystem/).
"""
struct CoupledSDEs{IIP, D, I, P} # <: ContinuousTimeDynamicalSystem
    integ::I
    # things we can't recover from `integ`
    p0::P
    diffeq # isn't parameterized because it is only used for display
    # D parametrised by length of u0
end

function CoupledSDEs(f, g, u0, p = SciMLBase.NullParameters();
        t0 = 0, diffeq = DEFAULT_DIFFEQ, noise_rate_prototype = nothing,
        noise = nothing, seed = UInt64(0))
    IIP = isinplace(f, 4) # from SciMLBase
    IIP == isinplace(g, 4) ||
        throw(ArgumentError("f and g must both be in-place or out-of-place"))

    s = correct_state(Val{IIP}(), u0)
    T = eltype(s)
    prob = SDEProblem{IIP}(f, g, s, (T(t0), T(Inf)), p,
        noise_rate_prototype = noise_rate_prototype, noise = noise, seed = seed)
    return CoupledSDEs(prob, diffeq)
end
# function CoupledSDEs(f, u0, p=SciMLBase.NullParameters(); σ, t0=0, diffeq=DEFAULT_DIFFEQ)
#     IIP = isinplace(f, 4) # from SciMLBase
#     CoupledSDEs(f, IIP ? idfunc! : idfunc, u0, p; t0=t0, diffeq=diffeq)
# end
# This preserves the referrenced MTK system and the originally passed diffeq kwargs
CoupledSDEs(ds::CoupledSDEs, diffeq) = CoupledSDEs(SDEProblem(ds), merge(ds.diffeq, diffeq))

function CoupledSDEs(prob::SDEProblem, diffeq = DEFAULT_DIFFEQ)
    IIP = isinplace(prob) # from SciMLBase
    D = length(prob.u0)
    P = typeof(prob.p)
    if prob.tspan === (nothing, nothing)
        # If the problem was made via MTK, it is possible to not have a default timespan.
        U = eltype(prob.u0)
        prob = SciMLBase.remake(prob; tspan = (U(0), U(Inf)))
    end
    solver, remaining = _decompose_into_solver_and_remaining(diffeq)
    integ = __init(prob, solver; remaining...,
        # Integrators are used exclusively iteratively. There is no reason to save anything.
        save_start = false, save_end = false, save_everystep = false,
        # DynamicalSystems.jl operates on integrators and `step!` exclusively,
        # so there is no reason to limit the maximum time evolution
        maxiters = Inf
    )
    return CoupledSDEs{IIP, D, typeof(integ), P}(integ, deepcopy(prob.p), diffeq)
end

# struct StochSystem
#     f::Function
#     g::Function
#     u::State
#     p::Parameters
#     σ::Float64
#     Σ::CovMatrix
#     process::Any

#     # Methods of StochSystem
#     function StochSystem(f, u, p, σ = 0.0, Σ = I(length(u)), process = "WhiteGauss")
#         new(f, idfunc, u, p, σ, Σ, process)
#     end
#     function StochSystem(f, g, u, p, σ = 0.0, Σ = I(length(u)), process = "WhiteGauss")
#         new(f, g, u, p, σ, Σ, process)
#     end
# end;

# Core functions acting on StochSystem
"""
    drift(sys::StochSystem, x::State)

Returns the drift field ``b(x)`` of the StochSystem `sys` at the state vector `x`.
"""
# drift(sys::StochSystem, x::State) = sys.f(x, sys.p, 0)

"""
    σg(sys::StochSystem)

Multiplies the noise strength `σ` of a StochSystem `sys` with its noise function `g`.
"""
function diag_noise_funtion(σ; in_place = false)
    if in_place
        return (du, u, p, t) -> σ .* idfunc!(du, u, p, t)
    else
        return (u, p, t) -> σ .* idfunc(u, p, t)
    end
end
function diag_noise_funtion(σ, g)
    if SciMLBase.isinplace(g, 4)
        return (du, u, p, t) -> σ .* g(du, u, p, t)
    else
        return (u, p, t) -> σ .* g(u, p, t)
    end
end

"""
    p(sys::StochSystem)
Concatenates the deterministic and stochastic parameter vectors `pf` and `pg` of a StochSystem `sys`.
"""
# function p(sys::StochSystem)
#     [sys.pf, sys.pg]
# end;

"""
    CoupledODEs(sys::StochSystem; diffeq, t0=0.0)

Converts a [`StochSystem`](@ref) into [`CoupledODEs`](https://juliadynamics.github.io/DynamicalSystems.jl/stable/tutorial/#DynamicalSystemsBase.CoupledODEs)
from DynamicalSystems.jl.
"""
# function CoupledODEs(
#         sys::StochSystem; diffeq = DynamicalSystemsBase.DEFAULT_DIFFEQ, t0 = 0.0)
#     DynamicalSystemsBase.CoupledODEs(
#         sys.f, SVector{length(sys.u)}(sys.u), sys.p; diffeq = diffeq, t0 = t0)
# end

# function to_cds(sys::StochSystem)
#     Base.depwarn("`to_cds` is deprecated, use `CoupledODEs(sys::StochSystem)` instead.",
#         :to_cds, force = true)
#     return CoupledODEs(sys)
# end

"""
    StochSystem(ds::CoupledODEs; g, σ, Σ, process)

Converts a [`CoupledODEs`](https://juliadynamics.github.io/DynamicalSystems.jl/stable/tutorial/#DynamicalSystemsBase.CoupledODEs)
system into a [`StochSystem`](@ref).
"""
# function StochSystem(ds::DynamicalSystemsBase.CoupledODEs; g = idfunc, σ = 0.0,
#         Σ = I(length(get_state(ds))), process = "WhiteGauss")
#     StochSystem(dynamic_rule(ds), g, current_state(ds), ds.p0, σ, Σ, process)
# end
