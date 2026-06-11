"""
    rate_track_return_tip(rs::RateSystem, Δts, Δps, mapper, ics; kw...)

Utilize the `global_continuation` functionality of Attractors.jl to
calculate a rate track-return-tip diagram for `rs` for a variety of
forcing duration and forcing scales `Δts, Δps` with the rate system starting always at `u0`.
Return:

1. a matrix of size `length(Δps)` × `length(Δts)` encoding the type of rate tipping behavior.
2. `attractors_cont`, the attractors of the global continuation of the unforced system.

## Keyword arguments

- `distance = Centroid()`: Distance function used when (1) matching attractors, and
  (2) mapping the end state of a rate simulation to its closest attractor through the
  `AttractorsViaProximity`.
- `proximity_kw`: Keywords propagated to [`AttractorsViaProximity`](@extref Attractors.AttractorsViaProximity), for mapping
  the rate system state to the unforced attractors at the middle and end (reverse) of forcing.
  Do not provide a `distance` as this is used from the above keyword.
- `decide_rate_outcome`: A three input function `f(track_id, return_id, start_id)`
  that returns an integer representing a possible rate tipping case of a simulation
  with a specific `Δt, Δp`. It inputs the attractor ID the system converges after the
  forwards rate forcing, the one it converges after the forcing is reverse, and the ID
  that `u0` converges at the starting parameter of the system.

## Description

This function formalizes and generalizes the concept of tracking, returning, or tipping,
in rate forced systems introduced in [Ritchie2023](@cite).
To achieve this, it uses global continuation using `mapper, ics`.
The function will run a whole `global_continuation` run over the
parameter range `prange = p0 .+ Δps` to establish the unfrozen system attractors and assign
unique IDs to them throughout `prange`.
If you instead want to run the global continuation yourself,
which allows you to change `Δts, Δps` without re-running the global continuation,
then simply do:
```julia
pcurve = unforced_pcurve(rs, Δps)
matcher = MatchBySSSetDistance(distance)
ascm = AttractorsSeedContinueMatch(mapper, matcher)
_, attractors_cont = global_continuation(ascm, pcurve, ics)
rate_track_return_tip(rs, Δts, Δps, attractors_cont; distance, kw...)
```

After the global continunation, the function will perform a
rate simulation with duration `Δt` and scale `Δp` for all
combinations.
For each, it will then assign an integer corresponding
to the type of rate-dependent behaviour as in [Ritchie2023](@cite).
By default, these integers are ∈ (1, 2, 3) and mean:

1. Tracking always.
2. Return but not tracking (safe overshoot).
3. Tipping always (failure to track).

You can however provide a custom function that may have more involved decision logic
via the keyword `decide_rate_outcome`.
To find the unforced system attractor IDs at the middle and and of the rate simulation
an `AttractorsViaProximity` is used.

## Notes

This function applies the same duration and scaling to all parameters that are forced in `rs`,
and enforces all forcings to be reversed as well.
The profiles of each parameter are individual though.
If any profile does not have the `reverse = true` option, an error is thrown.
"""
function rate_track_return_tip(
        rs::RateSystem, Δts, Δps, mapper::Attractors.AttractorMapper, ics;
        distance = StateSpaceSets.Centroid(), kw...
    )
    pcurve = unforced_pcurve(rs, Δps)
    matcher = Attractors.MatchBySSSetDistance(; distance)
    ascm = Attractors.AttractorSeedContinueMatch(mapper, matcher)
    _, attractors_cont = Attractors.global_continuation(ascm, pcurve, ics)
    return rate_track_return_tip(rs, Δts, Δps, attractors_cont; distance, kw...)
end

function rate_track_return_tip(
        rs::RateSystem, Δts, Δps, attractors_cont::AbstractVector;
        distance = StateSpaceSets.Centroid(),
        proximity_kw = NamedTuple(),
        u0 = initial_state(rs),
        decide_rate_outcome = decide_rate_outcome_default,
    )
    if any(isequal(false), values(rs.specs.forcing_reverse))
        error("Provided `RateSystem` must have a `reverse=true` option for all profiles.")
    end
    # configure proximity, finding starting attractor
    pcurve = unforced_pcurve(rs, Δps)
    unforced = rs.specs.unforced_system
    function find_attractor_id(unforced, u, p, attractors)
        set_parameters!(unforced, p)
        proximity = Attractors.AttractorsViaProximity(unforced, attractors; proximity_kw..., distance)
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
            # run forwards forcing. `step!` advances by a *relative* time, and `reinit!`
            # resets to the system's initial time, so step to the end of the forward
            # interval [tstart, tstart + dt] relative to that reset time.
            DynamicalSystemsBase.reinit!(rs, u0)
            t_reinit = current_time(rs)
            SciMLBase.step!(rs, (tstart + dt) - t_reinit, true)
            u = current_state(rs)
            track_id = find_attractor_id(unforced, u, pcurve[i], attractors_cont[i])
            # run backwards forcing over the reverse interval of equal duration
            SciMLBase.step!(rs, dt, true)
            u = current_state(rs)
            return_id = find_attractor_id(unforced, u, pcurve[1], attractors_cont[1])
            # deduce tipping type
            rate_type[i, j] = decide_rate_outcome(track_id, return_id, start_id)
        end
    end
    return rate_type, attractors_cont
end

function decide_rate_outcome_default(track_id, return_id, start_id)
    if track_id == return_id == start_id
        return 1 # track always
    elseif return_id == start_id && track_id ≠ start_id
        return 2 # return but not track
    elseif return_id ≠ start_id && track_id ≠ start_id
        return 3 # always tip
    else
        error("Unclassified rate outcome (track_id=$track_id, return_id=$return_id, \
               start_id=$start_id); provide a custom `decide_rate_outcome`.")
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
