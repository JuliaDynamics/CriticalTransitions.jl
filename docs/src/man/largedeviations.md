# Large deviation theory

## Action functionals
```@docs
fw_action(sys::StochSystem, path, time)
om_action(sys::StochSystem, path, time)
```

## Calculating instantons (minimum action paths)

```@docs
mam(sys::StochSystem, x_i::State, x_f::State, N::Int, T::Float; kwargs...)
```