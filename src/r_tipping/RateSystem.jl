# Consider the ODE dxₜ/dt = f(xₜ,p). We want to ramp one parameter of this ODE.

# We ask the user to define: 
#  1) an autonomous ContinuousTimeDynamicalSystem
#  2) a protocol for the time-dependent forcing via the struct RateConfig

# Then we give back the ContinuousTimeDynamicalSystem with the parameter 
# changing according to the rate protocol 'RateConfig'
"""
    RateConfig

Time-dependent forcing protocol ``p(t)`` describing the evolution of a parameter over a
time interval `[section_start, section_end]`.
Used to construct a non-autonomous `RateSystem`.

## Fields
- `pfunc`: function ``p(t)`` from ``R → R`` describing the parameter time dependence
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
    RateSystem(sys::ContinuousTimeDynamicalSystem, rc::RateConfig, pidx; kwargs...)

Creates a `RateSystem` type from an autonomous dynamical system `sys` and a time-dependent
parametric forcing protocol of `RateConfig` type.

## Keyword arguments
- `forcing_start = RateConfig.section_start`: Time when parameter shift starts (before this, the resulting system will be autonomous)
- `forcing_length = RateConfig.section_end - RateConfig.section_start`: Time-interval over which RateConfig.pfunc([RateConfig.section_start, RateConfig.section_end]) is spread out (for window_length > RateConfig.section_end - RateConfig.section_start) or squeezed into (for window_length < RateConfig.section_end - RateConfig.section_start)
- `forcing_scale = 1.0`: Amplitude of the ramping. The ramping is then automatically rescaled 
- `t0 = 0.0`: Initial time of the resulting non-autonomous system (relevant to later compute trajectories)

## Fields of the resulting `RateSystem` type
- `system`: Nonautonomous [`CoupledODEs`](@ref) constructed from an autonomous system `sys` and a `RateConfig`
- `forcing`: Function giving the value of the ramped parameter for each time.
- `pidx`: Index of the ramped parameter within the parameter container of the autonomous system

The returned `RateSystem.system` is is 
autonomous before `forcing_start`, 
non-autnonmous from `forcing_start` to `forcing_start+forcing_length` with the parameter shift given by the [`RateConfig`](@def), and again 
autonomous after `forcing_start+forcing_length`.
"""
function RateSystem(auto_sys::ContinuousTimeDynamicalSystem, rc::RateConfig, pidx;
    forcing_start = rc.section_start,
    forcing_length = rc.section_end - rc.section_start,
    forcing_scale = 1.0,
    t0 = 0.0)

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


