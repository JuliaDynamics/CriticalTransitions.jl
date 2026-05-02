# =====================================================================
# Time evolution of probability densities under the discretised
# Fokker–Planck equation.
# =====================================================================

"""
$(TYPEDSIGNATURES)

Time-evolved probability density. Solves the discrete Fokker–Planck
equation `dρ/dt = Qᵀ ρ` (where `Qᵀ = `[`fokker_planck_operator`](@ref)`(gen)`)
starting from `ρ_0`, and records `ρ` on a uniform time grid.

Returns `(ρs, t)` where `t = Ttr:Δt:(Ttr + T)` and `ρs` is a
`Matrix{Float64}` of size `(ncells(gen), length(t))` whose `i`-th column
is `ρ(t[i])`.

The signature mirrors `DynamicalSystemsBase.trajectory`:

```julia
ρs, t = propagate_density(gen, T, ρ_0)                      # snapshot at 0 and T
ρs, t = propagate_density(gen, T, ρ_0; Δt = 0.1)            # snapshots every 0.1
ρs, t = propagate_density(gen, T, ρ_0; Δt = 0.1, Ttr = 5)   # skip transient
```

# Arguments
- `gen` — diffusion generator (provides `Qᵀ`).
- `T` — total time over which the trajectory is recorded.
- `ρ_0` — initial density (length-`ncells(gen)` vector).

# Keyword arguments
- `Δt = T` — sampling interval. Default returns just the initial and
  final snapshots.
- `Ttr = 0` — transient time skipped before recording starts.
- `tol = 1e-7` — Krylov local-error tolerance.
- `m = 30` — maximum Krylov subspace dimension before restart.
- `adaptive = true` — use adaptive Krylov substepping. Disable only for
  debugging or to force fixed-step Krylov of dimension `m`.

For mass-preserving boundary conditions ([`Reflecting`](@ref) /
[`Periodic`](@ref)) the total mass is preserved across the trajectory;
for [`Absorbing`](@ref) boundaries it decays as probability leaks
through. Internal substepping is automatic — `tol`/`m`/`adaptive`
control the Krylov-subspace matrix-exponential algorithm in
[`ExponentialUtilities`](https://docs.sciml.ai/ExponentialUtilities/stable/),
which handles stiffness without ever forming `exp(t·Qᵀ)` explicitly.
"""
function propagate_density(
    gen::DiffusionGenerator,
    T::Real,
    ρ_0::AbstractVector;
    Δt::Real=T,
    Ttr::Real=0,
    tol::Real=1e-7,
    m::Integer=30,
    adaptive::Bool=true,
)
    length(ρ_0) == ncells(gen) || throw(
        DimensionMismatch(
            "ρ_0 has length $(length(ρ_0)) but generator has $(ncells(gen)) cells",
        ),
    )
    T >= 0 || throw(ArgumentError("T must be ≥ 0; got $T"))
    Δt > 0 || throw(ArgumentError("Δt must be > 0; got $Δt"))
    Ttr >= 0 || throw(ArgumentError("Ttr must be ≥ 0; got $Ttr"))

    t = collect(Float64, Ttr:Δt:(Ttr + T))
    F = fokker_planck_operator(gen)
    b = Vector{Float64}(ρ_0)

    # expv_timestep returns a Vector for a scalar / single-element t and a
    # Matrix for length(t) ≥ 2. Normalise to (N, length(t)) Matrix output.
    ρs = if length(t) == 1
        reshape(
            expv_timestep(t[1], F, b; tol=tol, m=m, adaptive=adaptive), :, 1,
        )
    else
        expv_timestep(t, F, b; tol=tol, m=m, adaptive=adaptive)
    end
    return ρs, t
end
