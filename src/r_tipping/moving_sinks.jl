"""
    moving_fixedpoints(rs::RateSystem, box; kwargs...)

Calculates fixed points of a nonautonomous dynamical system at snapshots in time, giving
insight into how the equilibria change as the forcing changes over time.

## Input
`rs` is a dynamical system of type [`RateSystem`](@ref)
and `box` is a vector of intervals (e.g. `[interval(0,1), interval(-1,1)]`)
that delimits the phase space volume in which to look for equilibria
(see [`Attractors.fixedpoints`](@ref)).

## Keyword arguments
- `Δt = 0.1`: Time step between fixed point calculations (in system time units)

## output
Returns a triple of vectors:
1. Fixed points found at each time point
2. Eigenvalues associated with the fixed points
3. Stability of the fixed points (1 - stable; 0 - unstable)

The term "moving sinks" refers to Wieczorek et al. (2023).
"""
function moving_fixedpoints(rs::RateSystem, box; Δt=0.1)
    fp, eig, stab = [], [], []
    t_start = rs.forcing.forcing_start_time
    t_end = rs.forcing_forcing_start_time + rs.forcing.forcing_duration
    times = t_start:Δt:t_end

    for t in times
        _fp, _eig, _stab = fixedpoints(frozen_system(rs, t), box)
        push!(fp, _fp)
        push!(eig, _eig)
        push!(stab, _stab)
    end
    return fp, eig, stab
end

# modified version
#function moving_sinks(
#    ds::ContinuousTimeDynamicalSystem, rp::RateProtocol, box; 
#    times=rp.t_start/rp.r:(rp.t_end/rp.r-rp.t_start/rp.r)/100:rp.t_end/rp.r
#)
#    fp, eig, stab = [], [], []
#    λ_mod  = rp.λ
#    t_start = rp.t_start
#    t_end = rp.t_end
#    for t in times
#        
#        t̃ = if r*t ≤ t_start
#                t_start
#            elseif t_start < r*t < t_end
#                r*t
#            else
#                t_end
#            end; # the value of the argument of the λ function in the drift field at time t
#        
#        # the fixed point function computes fixed points for the autonomous system attained at time t = 0
#        # ensuring that rp.t_start < 0 < rp.t_end and shifting the rp.λ function so that modified drift calls λ_mod(p,t̃) = λ_mod(p,r*0+t̃), rather than λ_mod(p,t_start+t̃) or λ_mod(p,t_end+t̃) for instance
#        rp.t_start = -1.0
#        rp.t_end = 1.0
#        rp.λ = (p,τ) -> λ_mod(p,τ+t̃)
#        
#        rate_sys = apply_ramping(ds, rp, rp.t_start)
#        _fp, _eig, _stab = fixedpoints(rate_sys, box)
#        push!(fp, _fp)
#        push!(eig, _eig)
#        push!(stab, _stab)
#    end
#    return fp, eig, stab
#end
