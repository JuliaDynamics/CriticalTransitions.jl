"""
    rate_track_return_tip(rs::RateSystem, Δts, Δps, mapper, u0; kw...)

Utilize the `global_continuation` functionality of Attractors.jl to
calculate a rate track-return-tip diagram for `rs` for a variety of
forcing duration and forcing scales `Δts, Δps` with the rate system starting always at `u0`.
Return:

1. a matrix of sides `Δts` × `Δps` encoding the type of rate tipping behavior.
2. `attractors_cont`, the attractors of the global continuation of the unforced system.

## Keyword arguments

- `distance = Centroid()`: Distance function used when (1) matching attractors, and
  (2) mapping the end state of a rate simulation to its closest attractor through the
  `AttractorsViaProximity`.
- `proximity_kw`: Keywords propagated to [`AttractorsViaProximity`](@ref), for mapping
  the rate system state to the unforced attractors at the middle and end (reverse) of forcing.
  At a minimum you'd want to pass here the same distance function as the one used during global
  continuation (which happens automatically for the default values).

## Description

This function formalizes and generalizes the concept of tracking, returning, or tipping,
in rate forced systems introduced in [Ritchie2023](@cite).
To achieve this, it uses global continuation, that must be run over the same parameter
range as the one covered by `Δps`:

```julia
pcurve = unforced_pcurve(ratestommel, Δps)

# you decide how to do global continuation:
mapper = some AttractorsMapper( ... )
matcher = MatchBySSSetDistance( some distance function ... )
sampler = some initial conditions sampler
ascm = AttractorsSeedContinueMatch(mapper, matcher)
_, attractors_cont = global_continuation(ascm, pcurve, sampler)
```
`attractors_cont`

The function will run a whole `global_continuation` run over the
parameter range `prange = p0 .+ Δps` to establish the unfrozen system attractors and assign
unique IDs to them throughout `prange`. If you have already done a global continuation,
use the call signature below.

After the continunation, it will perform a rate simulation with duration `Δt` and scale `Δp` for all
combinations, assigning to each combination an integer ∈ (1, 2, 3). These correspond
to the type of rate-dependent behaviour as in [Ritchie2023](@cite), as:

1. Tracking always.
2. Return but not tracking (safe overshoot).
3. Tipping always (failure to track).

To decide which of the three occurs, `AttractorsViaProximity` is used, mapping the state
to the attractor it converges and checking if that coincides with the starting one.

## Notes

This function applies the same duration and scaling to all parameters that are forced in `rs`,
and enforces all forcings to be reversed as well.
The profiles of each parameter are individual though.
"""
function rate_track_return_tip(
        rs::RateSystem, Δts, Δps, mapper::AttractorsMapper, ics;
        distance = Centroid(), kw...
    )
    pcurve = unforced_pcurve(rs, Δps)
    matcher = MatchBySSSetDistance(distance)
    ascm = AttractorsSeedContinueMatch(mapper, matcher)
    _, attractors_cont = global_continuation(ascm, pcurve, ics)
    return rate_track_return_tip(rs, Δts, Δps, attractors_cont; u0, distance, proximity_kw)
end

function rate_track_return_tip(
        rs::RateSystem, Δts, Δps, attractors_cont::AbstractVector;
        distance = Centroid(),
        proximity_kw = NamedTuple(),
        u0 = initial_state(rs),
    )
    if any(isfalse, values(rs.spec.forcing_reverse))
        error("Provided `RateSystem` must have a `reverse=true` option for all profiles.")
    end
    # configure proximity, finding starting attractor
    unforced = rs.specs.unforced_system
    function find_attractor_id(unforced, u, p, attractors)
        set_parameters!(unforced, p)
        proximity = AttractorsViaProximity(unforced, attractors; proximity_kw..., distance)
        return proximity(u)
    end
    start_id = find_attractor_id(unforced, u0, pcurve[1], attractors_cont[1])

    # find critical parameter (although we don't use this anymore):
    # critical parameter = when attractor doesn't exist anymore
    pc_index = findlast(d -> haskey(d, start_id), attractors_cont)
    pc = pcurve[pc_index]

    # homogenize start time
    set_forcing_start!(rs, first(values(rs.specs.forcing_start_time)))
    tstart = first(values(rs.specs.forcing_start_time))

    # run the main computation
    rate_type = zeros(Int, length(Δps), length(Δts))

    for (i, dp) in enumerate(Δps)
        set_forcing_scale!(rs, dp)
        for (j, dt) in enumerate(Δts)
            set_forcing_duration!(rs, dt)
            # run forwards forcing
            DynamicalSystemsBase.reinit!(rs, u0)
            T = tstart + dt
            step!(rs, T, true)
            u = current_state(rs)
            track_id = find_attractor_id(unforced, u, pcurve[i], attractors_cont[i])
            # run backwards forcing
            step!(rs, dt, true)
            u = current_state(rs)
            return_id = find_attractor_id(unforced, u, pcurve[1], attractors_cont[1])
            # deduce tipping type
            rate_type[i, j] = label_rate_outcome(track_id, return_id, start_id)
        end
    end
    return rate_type
end

function label_rate_outcome(track_id, return_id, start_id)
    if track_id == return_id == start_id
        return 1
    elseif return_id == start_id && track_id ≠ start_id
        return 2 # return but not track
    elseif return_id ≠ start_id && track_id ≠ start_id
        return 3 # always tip
    end
end


# TODO: Generalize this to Dps being a vector of dictionaries.
"""
    unforced_pcurve(rs::RateSystem, Δps::AbstractVector)

Given a range of parameter increments `Δps`, return a parameter curve `pcurve`
corresponding to all parameter values that `rs` will be forced through.
This `pcurve` can be used to perform a global continuation for the unforced system.
"""
function unforced_pcurve(rs::RateSystem, Δps)
    p0s = rs.specs.p0
    forced_pkeys = keys(rs.specs.forcing_profile)
    return [Dict(pidx => current_parameter(rs, pidx, p0s) + p for pidx in forced_pkeys) for p in Δps]
end


