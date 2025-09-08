# we consider the ODE dxₜ/dt = f(xₜ,p(rt))
# here p = p(t) ∈ Rᵐ is a function containing all the system parameters 

# We ask the user to define: 
#  1) a ContinuousTimeDynamicalSystem that should be investigated and
#  2) a protocol for the time-dependent forcing with the struct RateConfig

# Then we give back the ContinuousTimeDynamicalSystem with the parameter 
# changing according to the rate protocol
"""
    RateConfig

Time-dependent forcing protocol containing the information to apply a parameter shift to an autonomous system.

Fields
==============

- pidx: index of the parameter vector from the underlying autonomous system which is made to be time-dependent
- p: monotonic function from R \rightarrow R (which describes the time-dependent parametric forcing) 
- section_start: time from which the parameter ramping function is considered
- section_end: time until which the parameter ramping function is considered
# add more in the docs about the rationale behind the construction of the nonautonomous system
- window_start: the time at which the nonautonomous dynamics starts for the resulting system  
- window_length: the duration of time of nonautonomous dynamics within the resulting system [i.e. within (window_start,window_start+window_length)] 
- dp: the difference in the explicitly time-dependent parameter value attained across the ramping

Default values 
==============

- t_start = -window_length/2
- dp = 1
"""
mutable struct RateConfig
    pidx::Int64
    p::Function  
    section_start::Float64
    section_end::Float64  
    window_start::Float64
    window_length::Float64
    dp::Float64
end

## convenience functions

RateConfig(pidx::Int64,p::Function,section_start::Float64,section_end::Float64,window_length::Float64,dp::Float64) = RateConfig(pidx,p,section_start,section_end,-window_length/2,window_length,dp)
RateConfig(pidx::Int64,p::Function,section_start::Float64,section_end::Float64,window_start::Float64,window_length::Float64) = RateConfig(pidx,p,section_start,section_end,window_start,window_length,1.0)
RateConfig(pidx::Int64,p::Function,section_start::Float64,section_end::Float64,window_start::Float64,window_length::Float64) = RateConfig(pidx,p,section_start,section_end,-window_length/2,window_length,1.0)

## the following function creates a piecewise constant function in alignment with the entries of the RateConfig and the parameter value of the underlying autonomous system
function p_modified(p0::Float64,rc::RateConfig,t::Float64)

    # extracting the entries of rc
    pidx = rc.pidx
    p = rc.p
    section_start = rc.section_start
    section_end = rc.section_end
    window_start = rc.window_start 
    window_length = rc.window_length 
    dp = rc.dp

    # making the function piecewise constant with range [p0,p0+dp]    
    q = if t ≤ window_start
            return p0
        else
            if t < window_start+window_length
                # performing the time shift corresponding to stretching/squeezing 
                time_shift = ((section_end-section_start)/window_length)*(t-window_start)+section_start 
                return p0 + dp*(p(time_shift)-p(section_start))/(p(section_end)-p(section_start))
            else 
                return p0 + dp
            end
        end

    return q
end

"""
    apply_ramping(sys::CoupledODEs, rc::RateConfig, t0=0.0; kwargs...)

Applies a time-dependent [`RateConfig`](@def) to a given autonomous deterministic dynamical system
`sys`, turning it into a non-autonomous dynamical system. The returned [`CoupledODEs`](@ref)
has the explicit parameter time-dependence incorporated.

The returned [`CoupledODEs`](@ref) is autonomous from `t_0` to `t_start`, 
then non-autnonmous from `t_start` to `t_start+t_ramp_length` with the parameter shift given by the [`RateConfig`](@def),
and again autonomous from `t_start+t_ramp_length` to the end of the simulation:

`t_0`  autonomous    `t_start`  non-autonomous   `t_start+t_ramp_length`  autonomous   `∞`

Computing trajectories of the returned [`CoupledODEs`](@ref) can then be done in the same way as for any other [`CoupledODEs`](@ref).
"""
function apply_ramping(auto_sys::CoupledODEs, rc::RateConfig, t0=0.0; kwargs...)
    # we wish to return a continuous time dynamical system with modified drift field

    f = dynamic_rule(auto_sys)
    params = current_parameters(auto_sys)
    p0 = params[rc.pidx]

    function f_new(u, p, t)
        pvalue = p_modified(p0,rc,t)
        if p isa Union{AbstractArray, AbstractDict}
            setindex!(p, pvalue, rc.pidx)
        else
            setproperty!(p, rc.pidx, pvalue)
        end
        return f(u, p, t)
    end

    prob = remake(referrenced_sciml_prob(auto_sys); f=f_new, p=params, tspan=(t0, Inf))
    nonauto_sys = CoupledODEs(prob, auto_sys.diffeq)
    return nonauto_sys
end

## unresolved comments
## window_start initial comment: the parameter values of the past limit system are given by p0 + dp* p(p_parameteters,window_start)
## window_end initial comment: the parameter values of the future limit system are given by p(p_parameters,window_end) 
## - p_parameters: the vector of length mp giving parameters which are associated with p

