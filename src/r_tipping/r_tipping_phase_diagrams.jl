using CriticalTransitions, Attractors, CairoMakie

function stommel_f(x, p, t)
    T, S = x
    η1, η2, η3 = p
    q = abs(T - S)
    return SVector(η1 - T - q * T, η2 - η3 * S - q * S)
end

p = [2.6, 1, 0.3]
stommel = CoupledODEs(stommel_f, [0.3, 0.2], p)
profile = ForcingProfile(x -> cos(x)^2, (-π / 2, 0.0))
rs = RateSystem(stommel, profile, 1; forcing_start_time = 50.0)

N = 51
Δps = range(0.1, 1; length = N)
Δts = 2 .^ range(-4, -0.5; length = N)

u0 = [2.348128197247146, 2.455397131357698] # steady state at η0 = 2.6

function unforced_pcurve(rs::RateSystem, Δps)
    p0s = current_parameters(rs.specs.unforced_system)
    forced_pkeys = keys(rs.specs.forcing_profile)
    return [
        Dict(pidx => current_parameter(rs.specs.unforced_system, pidx, p0s) + p for pidx in forced_pkeys)
        for p in Δps
    ]
end

pcurve = unforced_pcurve(rs, Δps)

grid = (range(0, 10; length = 201), range(0, 10; length = 201))
ics, = statespace_sampler(grid)
mapper = AttractorsViaRecurrences(stommel, grid)
gca = AttractorSeedContinueMatch(mapper)
_, attractors_cont = global_continuation(
    gca, pcurve, ics; samples_per_parameter = 100
)

proximity_kw = (distance = Centroid(), ε = 0.1)

function attractor_id(unforced, u, p, attractors; proximity_kw)
    set_parameters!(unforced, p)
    proximity = AttractorsViaProximity(unforced, attractors; proximity_kw...)
    return proximity(u)
end

"""
    rate_track_return_tip(rs::RateSystem, Δts, Δps, attractors_cont, u0; kw...)

Utilize the `global_continuation` functionality of Attractors.jl to
calculate a rate track-return-tip diagram for `rs` for a variety of
forcing duration and forcing scales `Δts, Δps` for the rate system starting always at `u0`.
Return:

1. a matrix of sides `Δts` × `Δps` encoding the type of rate tipping behavior.
2. `attractors_cont`, the attractors global continuation of the frozen system.

## Keyword arguments

- `proximity_kw`: Keywords propagated to [`AttractorsViaProximity`](@ref), for mapping
  the rate system state to the unforced attractors at the middle and end (reverse) of forcing.

## Description

This function formalizes and generalizes the concept of tracking, returning, or tipping,
in rate forced systems introduced in [Ritchie2025](@cite).
To achieve this, it uses global continuation.

The function will run a whole `global_continuation` run over the
parameter range `prange = p0 .+ Δps` to establish the unfrozen system attractors and assign
unique IDs to them throughout `prange`. If you have already done a global continuation,
use the call signature below.

After the continunation, it will perform a rate simulation with duration `Δt` and scale `Δp` for all
combinations, assigning to each combination an integer ∈ (1, 2, 3). These correspond
to the type of rate-dependent behaviour as in [Ritchie2025](@cite), as:

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
    rs::RateSystem, Δts, Δps, attractors_cont, u0;
    proximity_kw = (distance = Centroid(),),
)
    unforced = rs.specs.unforced_system
    start_id = attractor_id(unforced, u0, pcurve[1], attractors_cont[1]; proximity_kw)

    set_forcing_start!(rs, first(values(rs.specs.forcing_start_time)))
    tstart = first(values(rs.specs.forcing_start_time))

    rate_type = zeros(Int, length(Δps), length(Δts))

    for (i, dp) in enumerate(Δps)
        set_forcing_scale!(rs, dp)
        for (j, r) in enumerate(Δts)
            dt = dp / r
            set_forcing_duration!(rs, dt)
            set_forcing_start!(rs, tstart)
            DynamicalSystemsBase.reinit!(rs.system, u0)
            tracktime = tstart + dt
            step!(rs, tracktime, true)
            u = current_state(rs)
            track_id = attractor_id(unforced, u, pcurve[i], attractors_cont[i]; proximity_kw)
            set_forcing_start!(rs, tracktime)
            set_forcing_scale!(rs, -dp)
            finaltime = tracktime + dt
            step!(rs, finaltime, true)
            u = current_state(rs)
            return_id = attractor_id(unforced, u, pcurve[1], attractors_cont[1]; proximity_kw)
            rate_type[i, j] = label_rate_outcome(track_id, return_id, start_id)
            set_forcing_scale!(rs, dp)
        end
    end
    return rate_type
end

function label_rate_outcome(track_id, return_id, start_id)
    if track_id == return_id == start_id
        return 1
    elseif return_id == start_id && track_id ≠ start_id
        return 2 # return but not track
    elseif return_id ≠ start_id
        return 3 # tip / fail to return
    end
end

rate_type = rate_track_return_tip(rs, Δts, Δps, attractors_cont, u0; proximity_kw)

start_id = attractor_id(rs.specs.unforced_system, u0, pcurve[1], attractors_cont[1]; proximity_kw)
Δp_bif = Δps[findlast(d -> haskey(d, start_id), attractors_cont)]

cmap = cgrad(["white", "#7c3aed", "#1e1b4b"], 3; categorical = true)
fig, ax, hm = heatmap(
    Δps, log2.(Δts), rate_type; colormap = cmap,
    colorrange = (0.5, 3.5)
)
ax.xlabel = "Δp"
ax.ylabel = "log(Δt)"
cb = Colorbar(fig[1, 2], hm)
cb.ticks = (1:3, ["always\ntrack", "return but\nnot track", "always\ntip"])
cb.ticklabelrotation = π / 2
fig
