# Define custom types
Parameters = Union{Vector{Any},Nothing};
CovMatrix = Union{Matrix,UniformScaling{Bool},Diagonal{Bool,Vector{Bool}}};
State = Union{Vector,SVector}

# Define RateSystem
"""
    RateSystem(f, pf, dim, σ, g, pg, Σ, process)
Defines a stochastic dynamical system with a rate dependent shift in `CriticalTransitions.jl`. See [documentation](https://juliadynamics.github.io/CriticalTransitions.jl/dev/man/stochsystem/).
"""
struct RateSystem
    # i think it would maybe be helpful to define a few more fields:
    # gL::Function, a function describing the noise functions of the time-dependent functions
    # dim2::Int64, gives the number of time-dependent functions
    # Σ₂::CovMatrix, gives the covariance matrix of the time-dependent parameters
    f::Function                 # vector of functions describing the derivatives of each state variable
    pf::Parameters              # parameters for the derivatives of the state variables
    td_inds::Vector{Bool}       # a boolean vector indicating which parameters are time-dependent
    L::Function                 # vector of functions describing the derivatives of each time-dependent parameter
    pL::Parameters              # parameters for the derivatives of the parameter variables
    T_trans::Float64            # the length of time before ramping begins
    T_shift::Float64            # the length of time where the specified parameters evolve
    dim::Int64
    σ::Float64
    g::Function
    pg::Parameters
    Σ::CovMatrix
    process::Any
end;

# Methods of RateSystem
function RateSystem(f, pf, L, T_trans, T_shift, dim)
    return RateSystem(
        f, pf, L, T_trans, T_shift, dim, 1.0, idfunc, nothing, I(dim), "WhiteGauss"
    )
end
function RateSystem(f, pf, L, T_trans, T_shift, dim, σ)
    return RateSystem(
        f, pf, L, T_trans, T_shift, dim, σ, idfunc, nothing, I(dim), "WhiteGauss"
    )
end
function RateSystem(f, pf, L, T_trans, T_shift, dim, σ, Σ)
    return RateSystem(f, pf, L, T_trans, T_shift, dim, σ, idfunc, nothing, Σ, "WhiteGauss")
end

# functions for RateSystem

"""
    simulate(sys::RateSystem, init; kwargs...)
Simulates the RateSystem `sys` forward in time, starting at initial condition `init`.
## Keyword arguments
* `dt=0.01`: time step of integration
* `tmax=1e3`: total time of simulation
* `solver=EM()`: numerical solver. Defaults to Euler-Mayurama
* `callback=nothing`: callback condition
* `progress=true`: shows a progress bar during simulation
* `kwargs...`: keyword arguments for `solve(SDEProblem)`
For more info, see [`SDEProblem`](https://diffeq.sciml.ai/stable/types/sde_types/#SciMLBase.SDEProblem).
> Warning: This function has only been tested for the `EM()` solver and out-of-place `SDEFunction`s.
"""
function simulate(
    sys::RateSystem,
    init;
    dt=0.01,
    tmax=1e3,
    solver=EM(),
    callback=nothing,
    progress=true,
    kwargs...,
)
    func(u, p, t) = fL(u, p, t, sys)

    prob = SDEProblem(func, σg(sys), init, (0, tmax), p(sys); noise=stochprocess(sys))

    return solve(prob, solver; dt=dt, callback=callback, progress=progress, kwargs...)
end;

"""
    relax(sys::StochSystem, init; kwargs...)
Simulates the deterministic dynamics of StochSystem `sys` in time, starting at initial condition `init`.

This function integrates `sys.f` forward in time, using the [`ODEProblem`](https://diffeq.sciml.ai/stable/types/ode_types/#SciMLBase.ODEProblem) functionality of `DifferentialEquations.jl`. Thus, `relax` is identical to [`simulate`](@ref) when setting the noise strength `sys.σ = 0`.

## Keyword arguments
* `dt=0.01`: time step of integration
* `tmax=1e3`: total time of simulation
* `solver=Euler()`: numerical solver. Defaults to explicit forward Euler
* `callback=nothing`: callback condition
* `kwargs...`: keyword arguments for `solve(ODEProblem)`

For more info, see [`ODEProblem`](https://diffeq.sciml.ai/stable/types/ode_types/#SciMLBase.ODEProblem).
For stochastic integration, see [`simulate`](@ref).

> Warning: This function has only been tested for the `Euler()` solver.
"""
function relax(
    sys::RateSystem,
    init;
    dt=0.01,
    tmax=1e3,
    solver=Euler(),
    callback=nothing,
    kwargs...,
)
    func(u, p, t) = fL(u, p, t, sys)

    prob = ODEProblem(func, init, (0, tmax), p(sys))
    return solve(prob, solver; dt=dt, callback=callback, kwargs...)
end;

"""
    σg(sys::RateSystem)
Multiplies the noise strength `σ` of a RateSystem `sys` with its noise function `g`.
"""
function σg(sys::RateSystem)
    # the parameter variation is purely deterministic (for now) - can implement noise functions for the time-dependent parameters later and change this function accordingly
    if is_iip(sys.f)
        function g_iip(du, u, p, t)
            return vcat(
                sys.σ .* sys.g(du[1:length(sys.u)], u[1:length(sys.u)], p, t),
                SVector{sum(sys.td_inds)}(zeros(sum(sys.td_inds))),
            )
        end
    else
        function g_oop(u, p, t)
            return vcat(
                sys.σ .* sys.g(u[1:length(sys.u)], p, t),
                SVector{sum(sys.td_inds)}(zeros(sum(sys.td_inds))),
            )
        end
    end
end;

"""
    p(sys::RateSystem)
Concatenates the deterministic and stochastic parameter vectors `pf`, 'pL', and `pg` of a RateSystem `sys`.
"""
function p(sys::RateSystem)
    return [vcat(sys.pf, sys.pL), sys.pg]
end;

"""
    stochprocess(sys::RateSystem)

Translates the stochastic process specified in `sys` into the language required by the
`SDEProblem` of `DifferentialEquations.jl`.
"""
function stochprocess(sys::RateSystem)
    if sys.process == "WhiteGauss"
        if sys.Σ == I(length(sys.u))
            return nothing
        else
            return gauss(sys)
        end
    else
        ArgumentError("Noise process not yet implemented.")
    end
end

"""
    gauss(sys::RateSystem)
Returns a Wiener process with dimension `length(sys.u)` and covariance matrix `sys.Σ`.

This function is based on the [`CorrelatedWienerProcess`](https://docs.sciml.ai/DiffEqNoiseProcess/stable/noise_processes/#DiffEqNoiseProcess.CorrelatedWienerProcess) of [`DiffEqNoiseProcess.jl`](https://docs.sciml.ai/DiffEqNoiseProcess/stable/), a component of `DifferentialEquations.jl`. The initial condition of the process is set to the zero vector at `t=0`.
"""
function gauss(sys::RateSystem)
    Σ₂ = zeros(sum(sys.td_inds), sum(sys.td_inds)) # the covariance matrix of the time-dependent parameters, fixed to zero for now
    # Returns a Wiener process for given covariance matrix and dimension of a StochSystem
    if is_iip(sys.f)
        W = CorrelatedWienerProcess!(
            [
                sys.Σ zeros(length(sys.u), sum(sys.td_inds))
                zeros(sum(sys.td_inds), length(sys.u)) Σ₂
            ],
            0.0,
            zeros(length(sys.u) + sum(sys.td_inds)),
        )
    else
        W = CorrelatedWienerProcess(
            [
                sys.Σ zeros(length(sys.u), sum(sys.td_inds))
                zeros(sum(sys.td_inds), length(sys.u)) Σ₂
            ],
            0.0,
            zeros(length(sys.u) + sum(sys.td_inds)),
        )
    end
    return W
end;

function fL(u, p, t, sys::RateSystem)
    pf = p[1][1:length(sys.pf)]
    pL = p[1][(length(sys.pf) + 1):end]

    #stationary = 0..t_trans ∪ t_trans+t_shift..Inf
    #du = t in stationary ? hcat(sys.f(u,p,t),zeros(length(sys.pf)) : hcat(sys.f(u,p,t),sys.L(u,p,t))

    # u is the vector made up of the state variables and time-dependent parameters
    state_vars = u[1:length(sys.u)] # the state-variables
    param_vars = u[(length(sys.u) + 1):end] # the (time-dependent) parameters

    pf[sys.td_inds] = param_vars # modifying the vector of fixed parameters to account for the current values of the time-dependent parameters

    nonstationary = sys.T_trans .. sys.T_trans + sys.T_shift # the time-interval of the ramping period

    return du = if t in nonstationary
        vcat(sys.f(state_vars, [pf], t), sys.L(param_vars, [pL], t))
    else
        vcat(sys.f(state_vars, [pf], t), SVector{sum(sys.td_inds)}(zeros(sum(sys.td_inds))))
    end
end;

function stochtorate(
    sys::StochSystem,
    td_inds::Vector{Bool},
    L::Function,
    pL::Vector,
    T_trans::Float64,
    T_shift::Float64,
)
    return RateSystem(
        sys.f,
        sys.pf,
        td_inds,
        L,
        pL,
        T_trans,
        T_shift,
        length(sys.u),
        sys.σ,
        sys.g,
        sys.pg,
        sys.Σ,
        sys.process,
    )
end
