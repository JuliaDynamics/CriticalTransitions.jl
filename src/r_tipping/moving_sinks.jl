"""
    moving_sinks(ds::ContinuousTimeDynamicalSystem, rp::RateProtocol, box; kwargs...)

Calculates fixed points of a nonautonomous dynamical system at snapshots in time, giving
insight into how the equilibria change as the forcing changes over time.

## Input
`ds` is a dynamical system of type `ContinuousTimeDynamicalSystem`, `rp` a `RateProtocol`,
and `box` is a vector of intervals (e.g. `[interval(0,1), interval(-1,1)]`)
that delimits the phase space volume in which to look for equilibria
(see [`fixedpoints`](@ref)).

## Keyword arguments
- `times = 0:0.1:1`: time points (relative to the system time `t`)

## output
Returns a triple of vectors:
1. Fixed points found at each time point
2. Eigenvalues associated with the fixed points
3. Stability of the fixed points (1 - stable; 0 - unstable)

The term "moving sinks" refers to Wieczorek et al. (2023).
"""
# function moving_sinks(
#     ds::ContinuousTimeDynamicalSystem, rp::RateProtocol, box; times=0:0.1:1
# )
#     fp, eig, stab = [], [], []
#     for t in times
#         rate_sys = apply_ramping(ds, rp, t)
#         _fp, _eig, _stab = fixedpoints(rate_sys, box)
#         push!(fp, _fp)
#         push!(eig, _eig)
#         push!(stab, _stab)
#     end
#     return fp, eig, stab
# end

# modified version
function moving_sinks(
    ds::ContinuousTimeDynamicalSystem, rp::RateProtocol, box; 
    times=rp.t_start/rp.r:(rp.t_end/rp.r-rp.t_start/rp.r)/100:rp.t_end/rp.r
)
    fp, eig, stab = [], [], []
    λ_mod  = rp.λ
    for t in times
        rp.λ = (p,τ) -> λ_mod(p,τ+t)
        rate_sys = apply_ramping(ds, rp, rp.t_start)
        _fp, _eig, _stab = fixedpoints(rate_sys, box)
        push!(fp, _fp)
        push!(eig, _eig)
        push!(stab, _stab)
    end
    return fp, eig, stab
end