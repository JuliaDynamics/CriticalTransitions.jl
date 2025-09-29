# we consider the ODE dxₜ/dt = f(xₜ,p(rt))
# here p = p(t) ∈ Rᵐ is a function containing all the system parameters 

# We ask the user to define: 
#  1) a ContinuousTimeDynamicalSystem that should be investigated and
#  2) a protocol for the time-dependent forcing via the struct RateConfig

# Then we give back the ContinuousTimeDynamicalSystem with the parameter 
# changing according to the rate protocol 'RateConfig'
"""
    RateConfig

Time-dependent forcing protocol ``p(t)`` describing the evolution of a parameter over a
time interval `[section_start, section_end]`.
Used to construct a non-autonomous `RateSystem`.

## Fields
- `pfunc`: function from ``R \rightarrow R`` describing the parameter time dependence
- `section_start`: start of the time interval over which `pfunc` is considered
- `section_end`: end of the time interval over which `pfunc` is considered

## Description
The `RateConfig` type allows to specify the functional form of a parametric
forcing over a time interval. This forcing protocol can then be applied to the parameter
of a dynamical system using the `RateSystem` constructor, which also allows to modify
the rate and magnitude of the forcing.
"""
mutable struct RateConfig
    pfunc::Function  
    section_start::Float64
    section_end::Float64  
end

"""
    RateSystem
"""
struct RateSystem
    system
    forcing
    pidx
end

"""
    RateSystem(sys::ContinuousDynamicalSystem, rc::RateConfig, pidx; kwargs...)

Creates a `RateSystem` type from an autonomous dynamical system `sys` and time-dependent
parametric forcing protocol of `RateConfig` type.
"""
function RateSystem(auto_sys::ContinuousTimeDynamicalSystem, rc::RateConfig, pidx;
    forcing_start=rc.section_start,
    forcing_length=rc.section_end - rc.section_start,
    forcing_scale=1.0,
    t0=0.0)

    params = deepcopy(current_parameters(auto_sys))
    p0 = params[pidx]

    if forcing_start < t0
        @warn "The forcing starts before the system initial time t0=$(t0)."
    end

    system = apply_ramping(   auto_sys, rc, pidx, p0, params, forcing_start, forcing_length, forcing_scale, t0)
    forcing = t -> p_modified(t, rc, p0, forcing_start, forcing_length, forcing_scale)

    return RateSystem(system, forcing, pidx)
end

## the following function creates a piecewise constant function in alignment with the entries of the RateConfig and the parameter value of the underlying autonomous system
function p_modified(t::Real, rc::RateConfig, p0, forcing_start, forcing_length, forcing_scale)

    p = rc.pfunc
    section_start = rc.section_start
    section_end = rc.section_end
    
    # making the function piecewise constant with range [p0,p0+forcing_scale]    
    q = if t ≤ forcing_start
            return p0
        else
            if t < forcing_start+forcing_length
                # performing the time shift corresponding to stretching/squeezing 
                time_shift = ((section_end-section_start)/forcing_length)*(t-forcing_start)+section_start 
                return p0 + forcing_scale*(p(time_shift)-p(section_start))/(p(section_end)-p(section_start))
            else 
                return p0 + forcing_scale
            end
        end

    return q
end

"""
    apply_ramping(sys::CoupledODEs, rc::RateConfig, t0=0.0; kwargs...)

Applies a time-dependent [`RateConfig`](@def) to a given autonomous deterministic dynamical system
`sys`, turning it into a non-autonomous dynamical system. The returned [`CoupledODEs`](@ref)
has the explicit parameter time-dependence incorporated.

The returned [`CoupledODEs`](@ref) is autonomous from `t_0` to `window_start`, 
then non-autnonmous from `window_start` to `window_start+window_length` with the parameter shift given by the [`RateConfig`](@def),
and again autonomous from `window_start+window_length` to the end of the simulation:

`t_0`  autonomous    `window_start`  non-autonomous   `window_start+window_length`  autonomous   `∞`

Computing trajectories of the returned [`CoupledODEs`](@ref) can then be done in the same way as for any other [`CoupledODEs`](@ref).
"""
function apply_ramping(auto_sys, rc, pidx, p0, params, forcing_start, forcing_length, forcing_scale, t0)
    # returns a continuous time dynamical system with modified drift field

    f = deepcopy(dynamic_rule(auto_sys))

    function f_new(u, p, t)
        pvalue = p_modified(t, rc, p0, forcing_start, forcing_length, forcing_scale)
        if p isa Union{AbstractArray, AbstractDict}
            setindex!(p, pvalue, pidx)
        else
            setproperty!(p, pidx, pvalue)
        end
        return f(u, p, t)
    end

    prob = remake(deepcopy(referrenced_sciml_prob(auto_sys)); f=f_new, p=params, tspan=(t0, Inf))
    nonauto_sys = CoupledODEs(prob, auto_sys.diffeq)
    return nonauto_sys
end


