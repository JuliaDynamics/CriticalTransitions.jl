using DynamicalSystemsBase: CoupledODEs, isinplace, __init, SciMLBase, correct_state, current_state, _set_parameter!
using StochasticDiffEq: SDEProblem, SDEIntegrator
using StochasticDiffEq: EM, SDEProblem

###########################################################################################
# DiffEq options
###########################################################################################
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


###########################################################################################
# Type
###########################################################################################

# Define CoupledSDEs
"""

"""
struct CoupledSDEs{IIP, D, I, P} <: ContinuousTimeDynamicalSystem
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

"""
    CoupledSDEs(ds::CoupledODEs; g, diffeq, noise_rate_prototype, noise, seed)

Converts a [`CoupledODEs`](https://juliadynamics.github.io/DynamicalSystems.jl/stable/tutorial/#DynamicalSystemsBase.CoupledODEs)
system into a [`CoupledSDEs`](@ref).
"""
function CoupledSDEs(ds::DynamicalSystemsBase.CoupledODEs; g,
        diffeq = DEFAULT_DIFFEQ, noise_rate_prototype = nothing,
        noise = nothing, seed = UInt64(0))
    CoupledSDEs(dynamic_rule(ds), g, current_state(ds), ds.p0; diffeq = diffeq,
        noise_rate_prototype = noise_rate_prototype, noise = noise, seed = seed)
end

"""
    CoupledODEs(sys::StochSystem; diffeq, t0=0.0)

Converts a [`CoupledSDEs`](@ref) into [`CoupledODEs`](https://juliadynamics.github.io/DynamicalSystems.jl/stable/tutorial/#DynamicalSystemsBase.CoupledODEs)
from DynamicalSystems.jl.
"""
function CoupledODEs(
        sys::CoupledSDEs; diffeq = DynamicalSystemsBase.DEFAULT_DIFFEQ, t0 = 0.0)
    DynamicalSystemsBase.CoupledODEs(
        sys.f, SVector{length(sys.u)}(sys.u), sys.p; diffeq = diffeq, t0 = t0)
end

# Pretty print
function additional_details(ds::CoupledSDEs)
    solver, remaining = _decompose_into_solver_and_remaining(ds.diffeq)
    return ["ODE solver" => string(nameof(typeof(solver))),
        "ODE kwargs" => remaining,
    ]
end

###########################################################################################
# API - obtaining information from the system
###########################################################################################

SciMLBase.isinplace(::CoupledSDEs{IIP}) where {IIP} = IIP
StateSpaceSets.dimension(::CoupledSDEs{IIP, D}) where {IIP, D} = D
DynamicalSystemsBase.current_state(ds::CoupledSDEs) = current_state(ds.integ)

function set_parameter!(ds::CoupledSDEs, args...)
    _set_parameter!(ds, args...)
    u_modified!(ds.integ, true)
    return
end

referrenced_sciml_prob(ds::CoupledSDEs) = ds.integ.sol.prob

# so that `ds` is printed
set_state!(ds::CoupledSDEs, u::AbstractArray) = (set_state!(ds.integ, u); ds)
SciMLBase.step!(ds::CoupledSDEs, args...) = (step!(ds.integ, args...); ds)

"""
    drift(sys::CoupledSDEs, x::State)

Returns the drift field ``b(x)`` of the CoupledSDEs `sys` at the state vector `x`.
"""
drift(sys::CoupledSDEs{IIP}, x) where {IIP} = IIP ? sys.f(x, sys.p0, 0) : sys.f(x, sys.p0, 0)