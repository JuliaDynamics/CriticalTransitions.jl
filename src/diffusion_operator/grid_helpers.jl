"""
$(TYPEDSIGNATURES)

Return a predicate for the open ball with given `center` and `radius`.
"""
ball(center, radius::Real) = x -> norm(x .- center) < radius

"""
$(TYPEDSIGNATURES)

Predicate matching points inside an axis-aligned box: returns the
closure `x -> all(lo .≤ x .≤ hi)`. `lo` and `hi` may be tuples,
`SVector`s, or any indexable types of the right length.
"""
cuboid(lo, hi) = x -> all(lo[k] <= x[k] <= hi[k] for k in eachindex(x))

"""
$(TYPEDSIGNATURES)

Predicate for the **sublevel set** `{x : f(x) < c}`. Useful for masking
by a potential, energy, or any scalar-valued function on state space.

```julia
using CriticalTransitions: sublevel
U(x) = 0.25 * (x[1]^2 - 1)^2 + 0.5 * x[2]^2
basin_low_energy = sublevel(U, 0.1)
```
"""
sublevel(f, c::Real) = x -> f(x) < c

"""
$(TYPEDSIGNATURES)

Reshape a per-cell vector (length `ncells(grid)`) back to a grid-shaped
`Array` of size `grid.nbox`. Useful for plotting and for spatial
operations on grid quantities.

```julia
using CriticalTransitions: reshape_to_grid
ρ = stationary_distribution(gen)
ρ_grid = reshape_to_grid(ρ, gen)
heatmap(grid.centers[1], grid.centers[2], ρ_grid)
```

Both [`DiffusionGenerator`](@ref) and [`CartesianGrid`](@ref) are
accepted as the second argument.
"""
function reshape_to_grid(v::AbstractVector, gen::DiffusionGenerator)
    return reshape_to_grid(v, gen.grid)
end
function reshape_to_grid(v::AbstractVector, grid::CartesianGrid)
    length(v) == ncells(grid) || throw(
        DimensionMismatch(
            "vector length $(length(v)) ≠ ncells = $(ncells(grid))",
        ),
    )
    return reshape(v, grid.nbox)
end
